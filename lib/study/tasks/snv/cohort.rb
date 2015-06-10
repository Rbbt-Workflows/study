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

  dep do |jobname,options|
    study = Study.setup(jobname.dup)
    jobs = study.genotyped_samples.collect{|sample| Sample.setup(sample, :cohort => study); sample.mutation_info(:job) }.flatten
    Misc.bootstrap(jobs, nil, :bar => "Processing sample mutation_info", :respawn => :always) do |job|
      job.produce
      nil
    end
    jobs
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

  dep do |jobname,options|
    study = Study.setup(jobname.dup)
    jobs = study.genotyped_samples.collect{|sample| Sample.setup(sample, :cohort => study); sample.gene_sample_mutation_status(:job) }.flatten
    Misc.bootstrap(jobs, nil, :bar => "Processing gene_sample_mutation_status", :respawn => :always) do |job|
      job.produce
      nil
    end
    jobs
  end
  task :sample_genes => :tsv do
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

  dep :sample_genes
  dep :genotyped_samples
  input :recurrent_threshold, :float, "Proportion of samples with gene affected (e.g. 0.05 for 5%)", 0.05
  task :recurrent_genes => :array do |recurrent_threshold|
    threshold_count = recurrent_threshold * step(:genotyped_samples).load.length
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

  dep :organism
  dep :genomic_mutations
  dep :num_genotyped_samples
  dep Sequence, :binomial_significance, :organism => :organism, :mutations => :genomic_mutations, :exome => true, :num_samples => :num_genotyped_samples
  task :binomial_significance => :tsv do
    Step.wait_for_jobs dependencies
    TSV.get_stream step(:binomial_significance)
  end

  dep :organism
  dep :genomic_mutations
  dep Sequence, :binomial_significance, :organism => :organism, :mutations => :genomic_mutations, :exome => :exome
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
    CMD.cmd("sort -k #{field_pos + 1} -g '#{step(:binomial_significance).path}' | cut -f 1", :pipe => true)
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
  dep Enrichment, :rank_enrichment, :organism => :organism, :list => :sorted_significant_genes, :background => :mappable_genes, :permutations => 100_000
  input :database, :string, "Database to use", nil, :select_options => Enrichment::DATABASES
  task :significance_rank_enrichment => :tsv do
    TSV.get_stream step(:rank_enrichment)
  end


  dep :mi
  dep :firestar
  dep :firestar_neighbours
  dep :interfaces
  dep :mutation_incidence
  dep :genomic_mutation_consequence
  dep :genotyped_samples
  task :firestar_analysis => :tsv do
    Step.wait_for_jobs dependencies
    ns_mis = Set.new step(:mi).load
    num_samples = step(:genotyped_samples).load.length
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
      total_mutations += samples.length
      mis.each do |mi|
        next unless firestar.include? mi or firestar_neighbours.include? mi or interfaces.include? mi
        protein, _sep, change = mi.partition(":")
        next unless Appris::PRINCIPAL_ISOFORMS.include? protein

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
    end

    set_info :total_mutations, total_mutations
    Open.write(file('broken_PPI'), affected_pairs.collect{|p| p*"\t"}*"\n")

    total_appris_sequence = self.total_appris_sequence
    feature_WGS_residues = self.feature_WGS_residues

    feature_info = {}
    feature_counts.each do |feature, sample_muts|
      count = sample_muts.uniq.length
      next if count < 3
      genes = sample_muts.collect{|sm| sm.split("_").last}.uniq
      wgs_count = feature_WGS_residues[feature]
      ratio = count.to_f / wgs_count.to_f
      feature_info[feature] = [ratio, count, wgs_count, genes*"|", sample_muts * "|"]
    end

    TSV.setup feature_info, :key_field => "Firestar feature", :fields => ["Ratio", "Matches", "Total residues with feature", "Ensembl Gene ID", "Sample Mutations"], :type => :double, :namespace => organism


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
      geom_text(vjust=0.5,hjust=-0.5) + stat_smooth(method="lm", se=FALSE) + 
      scale_x_log10() + scale_y_log10() + geom_point()
    EOF
    R::SVG.ggplotSVG(feature_info, script, 10, 10, :R_open => 'na.string=NULL')
  end

end
