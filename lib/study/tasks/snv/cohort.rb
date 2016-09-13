Workflow.require_workflow "MutationSignatures"
require 'study/tasks/mutation_signatures'

CohortTasks = Proc.new do

  helper :feature_WGS_residues do
    Persist.persist("WGS firestar feature counts", :marshal) do
      feature_WGS_residues = {}
      TSV.traverse Structure.appris_dataset do |transcript, values|
        transcript = transcript.first if Array === transcript
        next unless Appris::PRINCIPAL_TRANSCRIPTS.include? transcript
        gene, fire, spade, thump, crash = values
        next if seen.include? gene
        fire.each do |fpos|
          pos, _sep, f = fpos.partition ":"
          feature_WGS_residues[f] ||= 0
          feature_WGS_residues[f] += 1
        end
      end
      feature_WGS_residues
    end
  end

  helper :total_appris_sequence do
    Persist.persist("Total appris sequence length", :integer) do
      ensp2sequence = Organism.protein_sequence(organism).tsv :persist => true, :unnamed => true
      enst2ensp = Organism.transcripts(organism).index :persist => true, :target => "Ensembl Protein ID", :fields => ["Ensembl Transcript ID"], :unnamed => true
      total = 0
      TSV.traverse Structure.appris_dataset do |transcript, values|
        transcript = transcript.first if Array === transcript
        next unless Appris::PRINCIPAL_TRANSCRIPTS.include? transcript
        protein = enst2ensp[transcript]
        sequence = ensp2sequence[protein]
        next if sequence.nil?
        total += sequence.length
      end
      total
    end
  end

  dep Sample, :mutation_info, :compute => :bootstrap do |jobname,options|
    study = Study.setup(jobname.dup)
    study.genotyped_samples.collect{|sample| Sample.setup(sample, :cohort => study) unless Sample === sample; sample.mutation_info(:job, options) }.flatten
  end
  task :mutation_info => :tsv do
    Step.wait_for_jobs dependencies
    parser = TSV::Parser.new dependencies.first
    fields = parser.fields
    header = TSV.header_lines(parser.key_field, parser.fields, parser.options)

    io = Misc.open_pipe do |sin|
      sin.puts header

      TSV.traverse dependencies, :type => :array, :bar => "Joining mutation_info from samples" do |job|
        sample = job.clean_name.split(":").last
        TSV.traverse job.path, :type => :array do |line|
          next if line =~ /^#/
            sin.puts line
        end
      end
    end

    Misc.sort_stream io
  end

  dep Sample, :gene_sample_mutation_status, :compute => :bootstrap do |jobname,options|
    study = Study.setup(jobname.dup)
    study.genotyped_samples.collect{|sample| Sample.setup(sample, :cohort => study) unless Sample === sample; sample.gene_sample_mutation_status(:job, options) }.flatten
  end
  task :sample_gene_mutations => :tsv do
    Step.wait_for_jobs dependencies
    parser = TSV::Parser.new dependencies.first
    fields = parser.fields
    fields.unshift "Sample"
    header = TSV.header_lines(parser.key_field, parser.fields, parser.options.merge(:type => :double))

    io = Misc.open_pipe do |sin|
      sin.puts header

      TSV.traverse dependencies, :type => :array do |job|
        sample = job.clean_name.split(":").last
        TSV.traverse job, :type => :array do |line|
          next if line =~ /^#/
            gene,*rest = line.split("\t")
          parts = [gene, sample]
          parts.concat rest
          sin.puts parts * "\t"
        end
      end
    end

    TSV.collapse_stream io
  end

  dep Sample, :gene_cnv_status, :compute => :bootstrap do |jobname,options|
    study = Study.setup(jobname.dup)
    if study.has_cnv?
      study.cnv_samples.collect{|sample| Sample.setup(sample, :cohort => study) unless Sample === sample; sample.gene_cnv_status(:job, options) }.flatten
    else
      []
    end
  end
  task :sample_gene_cnvs => :tsv do
    if dependencies.any?
      parser = TSV::Parser.new dependencies.first
      fields = parser.fields
      fields.unshift "Sample"
      header = TSV.header_lines(parser.key_field, parser.fields, parser.options.merge(:type => :double))

      io = Misc.open_pipe do |sin|
        sin.puts header

        TSV.traverse dependencies, :type => :array do |job|
          sample = job.clean_name.split(":").last
          TSV.traverse job, :type => :array do |line|
            next if line =~ /^#/
              gene,*rest = line.split("\t")
            parts = [gene, sample]
            parts.concat rest
            sin.puts parts * "\t"
          end
        end
      end

      TSV.collapse_stream io
    else
      ""
    end
  end

  dep :sample_gene_mutations, :compute => :produce
  dep :sample_gene_cnvs, :compute => :produce
  task :sample_genes => :tsv do
    if study.has_cnv?
      io = TSV.paste_streams [step(:sample_gene_mutations), step(:sample_gene_cnvs)]
      parser = TSV::Parser.new io
      dumper = TSV::Dumper.new parser.options.merge(:fields => parser.fields[0..-3] + parser.fields[-1..-1])
      dumper.init
      TSV.traverse parser, :into => dumper do |gene,values|
        gene = gene.first if Array === gene

        samples, *rest = values
        cnv = rest.pop
        cnv_samples = rest.pop
        new_values = rest
        new_cnv = ['normal'] * samples.length
        cnv_samples.each_with_index do |cnv_sample,i|
          index = samples.index cnv_sample
          if index.nil?
            samples << cnv_sample
            new_values.each{|l| l << 'false'}
            new_cnv << cnv[i]
          else
            new_cnv[index] = cnv[i]
          end
        end
        new_values.unshift(samples)
        new_values.push(new_cnv)
        [gene, new_values]
      end
    else
      TSV.get_stream step(:sample_gene_mutations)
    end
  end

  dep :sample_genes
  input :recurrent_threshold, :float, "Proportion of samples with gene affected (e.g. 0.05 for 5%)", 0.05
  returns "Ensembl Gene ID"
  task :recurrent_genes => :array do |recurrent_threshold|
    threshold_count = (recurrent_threshold.to_f * study.genotyped_samples.length).ceil
    threshold_count = 2 if threshold_count < 2
    TSV.traverse step(:sample_genes), :type => :array, :into => :stream do |line|
      next if line =~ /^#/

      gene, samples, overlapping, affected = line.split("\t")
      affected_samples = samples.split("|").zip(affected.split("|")).collect{|d,a| a.to_s == 'true' ? d : nil}.compact
      next unless affected_samples.length >= threshold_count
      gene
    end
  end

  dep :mutation_incidence
  input :threshold, :integer, "Minimum number of samples with mutation", 2
  returns "Genomic Mutation"
  task :recurrent_mutations => :array do |threshold|
    TSV.traverse step(:mutation_incidence), :type => :flat, :into => :stream do |mut, samples|
      next unless samples.length >= threshold
      Array === mut ? mut.first : mut
    end
  end

  #dep :mutation_incidence
  #dep :genomic_mutation_consequence
  #task :mi_incidence => :tsv do |threshold|
  #  
  #  pasted = TSV.paste_streams [step(:genomic_mutation_consequence), step(:mutation_incidence)], :fix_flat => true

  #  dumper = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Sample"], :namespace => organism, :type => :flat
  #  dumper.init
  #  TSV.traverse pasted, :into => dumper, :type => :array, :bar => "Mutated Isoform incidence" do |line|
  #    next if line =~ /^#/
  #    mut, mis_str, sample_str = line.split("\t")

  #    res = []
  #    res.extend MultipleResult
  #    samples = sample_str.split("|") 
  #    mis_str.split("|").each do |mi|
  #      res << [mi, samples]
  #    end
  #    res
  #  end

  #  TSV.collapse_stream(dumper.stream)
  #end

  dep Sample, :mi, :file => nil, :vcf => false, :compute => :bootstrap do |jobname, options|
    study = Study.setup(jobname.dup)
    study.genotyped_samples.collect{|sample| Sample.setup(sample, :cohort => study) unless Sample === sample; sample.mi(:job, options) }.flatten
  end
  task :mi_incidence => :tsv do |threshold|
    
    FileUtils.mkdir files_dir unless File.exists? files_dir

    log :add_sample, "Add sample to mi result"
    all_file = file('all').find
    dependencies.each do |dep|
      io = TSV.get_stream dep
      sample = dep.clean_name.split(":").last
      tmp_file = file(sample).find
      Misc.consume_stream(io, false, tmp_file)
      `sed 's/$/\t#{sample}/' '#{tmp_file}' >> '#{all_file}'`
    end

    pasted = Misc.collapse_stream Misc.sort_stream(all_file)

    log :join, "Join pasted streams"
    dumper = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Sample"], :namespace => organism, :type => :double
    dumper.init
    TSV.traverse pasted, :into => dumper, :type => :array, :bar => "Mutated Isoform incidence" do |line|
      next if line =~ /^#/
      mi_str, sample_str = line.split("\t")
      [mi_str, sample_str.split("\t")]
    end

  end

  dep :organism
  dep :genomic_mutations, :compute => :produce
  dep :num_genotyped_samples
  dep Sequence, :binomial_significance, :organism => :organism, :mutations => :genomic_mutations, :exome => true, :num_samples => :num_genotyped_samples, :vcf => false
  task :binomial_significance => :tsv do
    Step.wait_for_jobs dependencies
    TSV.get_stream step(:binomial_significance)
  end

  dep :organism
  dep :genomic_mutations
  dep Sequence, :binomial_significance, :organism => :organism, :mutations => :genomic_mutations, :exome => :exome, :vcf => false
  task :binomial_significance_all => :tsv do
    Step.wait_for_jobs dependencies
    TSV.get_stream step(:binomial_significance)
  end

  dep :binomial_significance
  input :threshold, :float, "P-value threshold", 0.05
  task :significant_genes => :array do |threshold|
    TSV.traverse step(:binomial_significance), :into => :stream do |gene,values|
      next unless values[-1].to_f <= threshold
      gene
    end
  end

  dep :binomial_significance
  task :sorted_significant_genes => :array do |threshold|
    Step.wait_for_jobs dependencies
    field_pos = TSV.parse_header(step(:binomial_significance).path).all_fields.index "p.value"
    CMD.cmd("env LC_ALL=C sort -k #{field_pos + 1} -g -  | cut -f 1 |grep -v '#'", :pipe => true, :in => Open.open(step(:binomial_significance).join.path, :nocache => true))
  end

  dep :gene_damage_bias
  input :threshold, :float, "P-value threshold", 0.05
  task :damage_biased_genes => :array do |threshold|
    TSV.traverse step(:gene_damage_bias), :into => :stream do |gene,values|
      next unless values[-1].to_f <= threshold
      gene
    end
  end

  dep :gene_damage_bias
  task :sorted_damage_biased_genes => :array do |threshold|
    Step.wait_for_jobs dependencies
    field_pos = TSV.parse_header(step(:gene_damage_bias).path).all_fields.index "p.value"
    CMD.cmd("env LC_ALL=C sort -k #{field_pos + 1} -g -  | cut -f 1 |grep -v '#'", :pipe => true, :in => Open.open(step(:gene_damage_bias).join.path, :nocache => true))
  end

  task :mappable_genes => :array do |threshold|
    mappable_regions_file = Study.find_study(study).mappable_regions
    if mappable_regions_file.exists?
      Sequence.job(:genes_at_ranges, clean_name, :ranges => mappable_regions_file).run.values.flatten.compact.uniq
    else
      []
    end
  end

  Workflow.require_workflow "MutationEnrichment"
  dep :organism
  dep :mutation_incidence
  dep :mappable_genes
  dep MutationEnrichment, :sample_pathway_enrichment, :organism => :organism, :mutations => :mutation_incidence, :background => :mappable_genes
  input :database, :string, "Database to use", nil, :select_options => MutationEnrichment::DATABASES
  task :sample_enrichment => :tsv do
    TSV.get_stream step(:sample_pathway_enrichment)
  end

  dep :organism
  dep :genomic_mutations
  dep :mappable_genes
  dep MutationEnrichment, :mutation_pathway_enrichment, :organism => :organism, :mutations => :genomic_mutations, :background => :mappable_genes
  input :database, :string, "Database to use", nil, :select_options => MutationEnrichment::DATABASES
  task :mutation_enrichment => :tsv do
    TSV.get_stream step(:mutation_pathway_enrichment)
  end

  Workflow.require_workflow "Enrichment"
  dep :organism
  dep :recurrent_genes
  dep :mappable_genes
  dep Enrichment, :enrichment, :organism => :organism, :list => :recurrent_genes, :background => :mappable_genes
  input :database, :string, "Database to use", nil, :select_options => Enrichment::DATABASES
  task :recurrent_gene_enrichment => :tsv do
    TSV.get_stream step(:enrichment)
  end

  dep :organism
  dep :significant_genes
  dep :mappable_genes
  dep Enrichment, :enrichment, :organism => :organism, :list => :significant_genes, :background => :mappable_genes
  input :database, :string, "Database to use", nil, :select_options => Enrichment::DATABASES
  task :significant_gene_enrichment => :tsv do
    TSV.get_stream step(:enrichment)
  end

  dep :organism
  dep :sorted_significant_genes
  dep :mappable_genes
  dep Enrichment, :rank_enrichment, :organism => :organism, :list => :sorted_significant_genes, :background => :mappable_genes, :permutations => 50_000
  input :database, :string, "Database to use", nil, :select_options => Enrichment::DATABASES
  task :significance_rank_enrichment => :tsv do
    TSV.get_stream step(:rank_enrichment)
  end

  dep :organism
  dep :damage_biased_genes
  dep :mappable_genes
  dep Enrichment, :enrichment, :organism => :organism, :list => :damage_biased_genes, :background => :mappable_genes
  input :database, :string, "Database to use", nil, :select_options => Enrichment::DATABASES
  task :damage_biased_gene_enrichment => :tsv do
    TSV.get_stream step(:enrichment)
  end

  dep :organism
  dep :sorted_damage_biased_genes
  dep :mappable_genes
  dep Enrichment, :rank_enrichment, :organism => :organism, :list => :sorted_damage_biased_genes, :background => :mappable_genes, :permutations => 50_000
  input :database, :string, "Database to use", nil, :select_options => Enrichment::DATABASES
  task :damage_bias_rank_enrichment => :tsv do
    TSV.get_stream step(:rank_enrichment)
  end



  dep :mi
  dep :firestar
  dep :firestar_neighbours
  dep :interfaces
  dep :mutation_incidence
  dep :genomic_mutation_consequence
  task :firestar_analysis => :tsv do
    Step.wait_for_jobs dependencies
    ns_mis = Set.new step(:mi).load
    num_samples = study.genotyped_samples.length
    firestar = step(:firestar).path.tsv :fields => ["Firestar site"], :type => :flat
    firestar_neighbours = step(:firestar_neighbours).path.tsv :fields => ["Firestar neighbour site"], :type => :flat
    interfaces = step(:interfaces).path.tsv :fields => ["Partner Ensembl Protein ID"], :type => :flat
    mutation_index = step(:mutation_incidence).path.tsv :persist => true, :persist_file => file('mutation_index')

    ensp2ensg = Organism.transcripts(organism).index :persist => true, :target => "Ensembl Gene ID", :fields => ["Ensembl Protein ID"], :unnamed => true

    affected_pairs = []
    feature_counts = {}
    nfeature_counts = {}
    total_mutations = 0
    total_size = 0
    TSV.traverse step(:genomic_mutation_consequence) do |mut,mis|
      samples = mutation_index[mut]
      mis = ns_mis & mis
      next if mis.empty?
      found = false
      mis.each do |mi|
        protein, _sep, change = mi.partition(":")
        next unless Appris::PRINCIPAL_ISOFORMS.include? protein
        found = true
        next unless firestar.include? mi or firestar_neighbours.include? mi or interfaces.include? mi

        fire = firestar[mi] || []
        nfire = firestar_neighbours[mi] || []
        partner = interfaces[mi] || []
        fire.uniq!
        nfire.uniq!
        partner.uniq!

        gene = ensp2ensg[protein]
        partner.each do |pprotein|
          pgene = ensp2ensg[pprotein]
          affected_pairs << [gene,pgene]
        end

        sample_mutations = samples.collect{|s| [s,mut,gene] * "_" }
        fire.each do |f|
          feature_counts[f] ||= []
          feature_counts[f].concat sample_mutations
        end
        nfire.each do |f|
          nfeature_counts[f] ||= []
          nfeature_counts[f].concat sample_mutations
        end

      end
      total_mutations += samples.length if found
    end

    total_appris_sequence = self.total_appris_sequence

    set_info :total_mutations, total_mutations
    set_info :total_appris_sequence, total_appris_sequence

    Open.write(file('broken_PPI'), affected_pairs.collect{|p| p*"\t"}*"\n")

    feature_WGS_residues = self.feature_WGS_residues

    feature_info = {}
    feature_counts.each do |feature, sample_muts|
      count = sample_muts.uniq.length
      next if count < 3
      genes = sample_muts.collect{|sm| sm.split("_").last}.uniq
      wgs_count = feature_WGS_residues[feature]
      ratio = count.to_f / wgs_count.to_f
      freq = wgs_count.to_f/total_appris_sequence
      feature_info[feature] = [ratio, freq, count, wgs_count, genes*"|", sample_muts * "|"]
    end

    TSV.setup feature_info, :key_field => "Firestar feature", :fields => ["Ratio", "Frequency of residue", "Matches", "Total residues with feature", "Ensembl Gene ID", "Sample Mutations"], :type => :double, :namespace => organism


    feature_info = feature_info.R <<-EOF, nil, :R_open => 'na.string=NULL'

data$Feature = rownames(data)
names(data) <- make.names(names(data))
data$Matches = as.numeric(data$Matches)
data$Total.residues.with.feature = as.numeric(data$Total.residues.with.feature)

total_mutations = #{R.ruby2R total_mutations}
total_sequence = #{R.ruby2R total_appris_sequence}  #* #{num_samples}

data$p.value = apply(data,1,function(v){ 
  matches = as.numeric(v["Matches"])
  total = as.numeric(v["Total.residues.with.feature"])
  freq = total/total_sequence
  binom.test(matches, total_mutations, freq, alternative='greater')$p.value
})
data
    EOF

    feature_info.type = :double
    feature_info.key_field = "Firestar feature"
    feature_info.fields = feature_info.fields.collect{|f| f.gsub('p.value', 'p-value').gsub('.',' ') }
    Open.write(file('feature_info'), feature_info.to_s)

    nfeature_info = {}


    nfeature_counts.each do |feature, sample_muts|
      count = sample_muts.uniq.length
      next if count < 3
      genes = sample_muts.collect{|sm| sm.split("_").last}.uniq
      wgs_count = feature_WGS_residues[feature]
      ratio = count.to_f / wgs_count.to_f
      nfeature_info[feature] = [ratio, count, wgs_count,genes*"|", sample_muts * "|"]
    end

    TSV.setup nfeature_info, :key_field => "Firestar neighbour feature", :fields => ["Ratio", "Matches", "Total residues with feature", "Ensembl Gene ID", "Sample Mutations"], :type => :double, :namespace => organism

    "Done"
  end

  dep :firestar_analysis
  extension :svg
  task :firestar_plot => :text do
    feature_info = step(:firestar_analysis).file(:feature_info).tsv

    script =<<-EOF

    names(data) <- make.names(names(data))
    data$Ratio = as.numeric(data$Ratio)
    data$p.value = as.numeric(data$p.value)
    data$Feature = rownames(data); 
    ggplot(aes(x=Matches,y=Total.residues.with.feature, colour=log10(Ratio), label=Feature, size=-log10(p.value)), data=data) + 
      geom_text(vjust=0.5,hjust=-2) + stat_smooth(method="lm", se=FALSE) + 
      scale_x_log10() + scale_y_log10() 
    EOF
    R::SVG.ggplotSVG(feature_info, script, 10, 10, :R_open => 'na.string=NULL')
  end

  dep :broken_ppi
  dep :mutation_incidence
  dep :genomic_mutation_consequence
  task :broken_ppi_incidence => :tsv do
    Step.wait_for_jobs dependencies
    consequence = step(:genomic_mutation_consequence).path.index :target => "Genomic Mutation"
    incidence = step(:mutation_incidence).load
    consequence.unnamed = true
    incidence.unnamed = true
    dumper = TSV::Dumper.new :key_field => "Sample", :fields => ["PPI Pair"], :type => :double, :namespace => organism
    dumper.init
    TSV.traverse step(:broken_ppi), :into => dumper do |mi,values|
      mi = mi.first if Array === mi
      mutations = consequence[mi]
      source,target = values
      samples = incidence.values_at *mutations
      res = samples.flatten.collect{|sample|
        [sample, [[source, target] * "~"]]
      }
      res.extend MultipleResult
      res
    end

    Misc.collapse_stream dumper.stream
  end

  dep :firestar
  dep :genomic_mutation_consequence
  dep :mutation_incidence
  task :firestar_incidence => :tsv do
    Step.wait_for_jobs dependencies
    consequence = step(:genomic_mutation_consequence).path.index :target => "Genomic Mutation"
    incidence = step(:mutation_incidence).load
    consequence.unnamed = true
    incidence.unnamed = true
    dumper = TSV::Dumper.new :key_field => "Sample", :fields => ["Firestar feature"], :type => :double, :namespace => organism
    dumper.init
    TSV.traverse step(:firestar), :into => dumper do |mi,values|
      mi = mi.first if Array === mi
      mutations = consequence[mi]
      feature, range = values
      samples = incidence.values_at *mutations
      protein = mi.partition(":").first
      res = samples.flatten.collect{|sample|
        [sample, [[protein, feature] * ":"]]
      }
      res.extend MultipleResult
      res
    end

    Misc.collapse_stream dumper.stream
  end

  dep :kinmut
  dep :genomic_mutation_consequence
  dep :mutation_incidence
  task :kinmut_incidence => :tsv do
    consequence = step(:genomic_mutation_consequence).path.index :target => "Genomic Mutation"
    incidence = step(:mutation_incidence).load
    consequence.unnamed = true
    incidence.unnamed = true
    dumper = TSV::Dumper.new :key_field => "Sample", :fields => ["Firestar feature"], :type => :double, :namespace => organism
    dumper.init
    TSV.traverse step(:kinmut), :into => dumper do |mi,values|
      mi = mi.first if Array === mi
      mutations = consequence[mi]
      prediction = values.first
      samples = incidence.values_at *mutations
      protein = mi.partition(":").first
      res = samples.flatten.collect{|sample|
        [sample, [[protein, prediction] * ":"]]
      }
      res.extend MultipleResult
      res
    end

    Misc.collapse_stream dumper.stream
  end

  dep Study, :genomic_mutation_consequence do |jobname, options|
    Study.job(:genomic_mutation_consequence, jobname, {:principal => true}.merge(options))
  end
  dep :mutation_incidence
  task :kinase_sample_mutations => :tsv do 
    Workflow.require_workflow "KinMut2"
    require 'kinase'
    kinases_uni = KinMut2.all_kinases
    kinases_ensembl = Organism.identifiers(organism).tsv(:key_field => "UniProt/SwissProt Accession", :fields => ["Ensembl Protein ID"], :type => :flat, :persist => true).chunked_values_at(kinases_uni).flatten
    incidence = step(:mutation_incidence).load
    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Mutated Isoform", "Sample"], :namespace => organism, :type =>:double
    dumper.init
    TSV.traverse step(:genomic_mutation_consequence), :into => dumper do |mut,mis|
      next if mis.nil?
      kmis = mis.select{|mi| kinases_ensembl.include? mi.partition(":").first}
      next if kmis.empty?
      samples = incidence[mut]

      [mut, [mis, samples]]
    end
  end

  dep :mutation_incidence, :compute => :produce
  dep :genomic_mutations, :compute => :produce
  dep :sequence_ontology, :principal => true, :compute => :produce
  dep :organism
  dep :annotate_DbSNP
  dep :mi_damaged
  dep Sequence, :genes, :positions => :genomic_mutations, :organism => :organism
  dep Sequence, :reference, :positions => :genomic_mutations, :organism => :organism
  dep Sequence, :type, :mutations => :genomic_mutations, :organism => :organism
  task :mutation_extended_summary => :tsv do
    Step.wait_for_jobs dependencies

    damaged = step(:mi_damaged).load
    organism = step(:organism).load

    all_fields = TSV.parse_header(step(:annotate_DbSNP)).all_fields
    fields = ["CAF", "RS ID"]
    pos = fields.collect{|f| all_fields.index(f) + 1 }
    pos.unshift 1
    pos.unshift 3

    io = TSV.paste_streams [step(:mutation_incidence), step(:reference), step(:type), step(:genes), step(:sequence_ontology),CMD.cmd("cut -f #{pos * ","} '#{step(:annotate_DbSNP).path}'", :pipe => true)] , :fix_flat => true

    name_index = Translation.index(organism, "Associated Gene Name", "Ensembl Gene ID")

    io2 = TSV.traverse io, :type => :array, :into => :stream, :bar => true do |line|
      begin
        if line =~ /^#/
          if line =~ /^#:/
            if line.include? ":namespace"
              line
            else
              line.sub("#:", "#: :namespace=#{organism}")
            end
          else
            mut, sample, type, ref, genes, mi, mi_so, mut_so, so, rsid, orig, caf  = line.split("\t")
            [mut, sample, type, ref, genes, "Associated Gene Name", mi, "Damaged", mi_so, mut_so, so, rsid, "Original Mutation", caf, "Mutated Allele Frequency"] * "\t"
          end
        else
          mut, sample, type, ref, genes, mi, mi_so, mut_so, so, rsid, orig, caf  = line.split("\t")

          if caf and not caf.empty?
            alt = mut.split(":")[2]
            oalt = orig.split(":")[2].split(",")
            apos = oalt.index(alt)

            maf = apos.nil? ? "-" : caf.split(",")[apos+1]

          else
            maf = ""
          end

          if mi and not mi.empty?
            damaged = mi.split("|").collect{|m| damaged.include?(mi).to_s} * "|"
          else
            damaged = ""
          end

          gene_names = name_index.values_at(*genes.split("|")) * "|"
          [mut, sample, type, ref, genes, gene_names, mi, damaged, mi_so, mut_so, so, rsid, orig, caf, maf] * "\t"
        end
      rescue
        Log.exception $!
        raise $!
      end
    end

    io2
  end

  dep :mi_incidence
  task :mi => :array do
    TSV.traverse step(:mi_incidence), :type => :flat, :into => :stream do |mi, samples|
      (mi * samples.length) * "\n"
    end
  end

  dep :mi_incidence, :compute => :produce
  dep :firestar, :compute => :produce
  task :sample_mi_firestar => :tsv do
    Step.wait_for_jobs dependencies
    TSV.paste_streams dependencies, :fix_flat => true, :all_match => true
  end

  dep :mi_incidence, :compute => :produce
  dep :interfaces, :compute => :produce
  task :sample_mi_interfaces => :tsv do
    TSV.paste_streams dependencies, :fix_flat => true, :all_match => true
  end

  dep :mi_incidence, :compute => :produce
  dep :kinmut, :compute => :produce
  task :sample_mi_kinmut => :tsv do
    TSV.paste_streams dependencies, :fix_flat => true, :all_match => true
  end

  dep Sample, :num_genomic_mutations, :compute => :bootstrap do |jobname,options|
    study = Study.setup(jobname.dup)
    study.genotyped_samples.collect do |sample|
      sample.num_genomic_mutations(:job)
    end
  end
  task :num_genomic_mutations => :integer do
    dependencies.inject(0){|acc,dep| acc += dep.load}
  end

end
