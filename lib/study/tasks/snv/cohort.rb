CohortTasks = Proc.new do

  dep do |jobname,options|
    study = Study.setup(jobname.dup)
    jobs = study.genotyped_samples.collect{|sample| Sample.setup(sample, :cohort => study); sample.gene_sample_mutation_status(:job) }.flatten
    Misc.bootstrap(jobs, nil, :bar => "Processing sample gene_mutation_status") do |job|
      job.clean if job.dirty?
      job.run(false) unless job.done? or job.started? 
      Step.wait_for_jobs job unless job.done?
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

  dep :mutation_incidence
  input :threshold, :integer, "Minimum number of samples with mutation", 2
  returns "Genomic Mutation"
  task :recurrent_mutations => :array do |threshold|
    TSV.traverse step(:mutation_incidence), :type => :flat, :into => :stream do |mut, samples|
      next unless samples.length >= threshold
      Array === mut ? mut.first : mut
    end
  end

  dep do |jobname,options|
    study = Study.setup(jobname.dup)
    jobs = study.genotyped_samples.collect{|sample| Sample.setup(sample, :cohort => study); sample.mutation_info(:job) }.flatten
    Misc.bootstrap(jobs, nil, :bar => "Processing sample mutation_info") do |job|
      job.clean if job.dirty?
      job.run(true) unless (job.done? or job.started?) 
      job.join unless job.done?
    end
    jobs
  end
  task :mutation_info => :tsv do
    Step.wait_for_jobs dependencies
    parser = TSV::Parser.new dependencies.first
    fields = parser.fields
    header = TSV.header_lines(parser.key_field, parser.fields, parser.options)

    io = Misc.open_pipe do |sin|

      TSV.traverse dependencies, :type => :array do |job|
        sample = job.clean_name.split(":").last
        TSV.traverse job, :type => :array do |line|
          next if line =~ /^#/
            sin.puts line
        end
      end
    end

    Misc.sort_stream io
  end

  dep :organism
  dep :exome
  dep :genomic_mutations
  dep Sequence, :binomial_significance, :organism => :organism, :mutations => :genomic_mutations, :exome => :exome
  task :binomial_significance => :tsv do
    Step.wait_for_jobs dependencies
    TSV.get_stream step(:binomial_significance)
  end

  Workflow.require_workflow "MutationEnrichment"
  dep :organism
  dep :mutation_incidence
  dep MutationEnrichment, :sample_pathway_enrichment, :organism => :organism, :mutations => :mutation_incidence
  input :database, :string, "Database to use", nil, :select_options => MutationEnrichment::DATABASES
  task :sample_enrichment => :tsv do
    TSV.get_stream step(:sample_pathway_enrichment)
  end

  #dep do |jobname,options|
  #  study = Study.setup(jobname.dup)
  #  jobs = study.genotyped_samples.collect{|sample| Sample.setup(sample, :cohort => study); sample.mutation_genes(:job) }.flatten
  #  Misc.bootstrap(jobs) do |job|
  #    job.run(true) unless (job.done? or job.started?) and not job.dirty?
  #    job.join unless job.done?
  #  end
  #  jobs
  #end
  #task :mutation_basic_info => :tsv do

  #  parser = TSV::Parser.new dependencies.first
  #  fields = parser.fields
  #  good_fields = fields - ['homozygous']
  #  good_field_pos = good_fields.collect{|f| fields.index(f) + 2 }

  #  header = TSV.header_lines(parser.key_field, parser.fields, parser.options)
  #  io = Misc.open_pipe do |sin|
  #    sin.puts header

  #    TSV.traverse dependencies, :type => :array do |job|
  #      TSV.traverse job, :type => :array do |line|
  #        next if line =~ /^#/
  #        sin.puts line
  #      end
  #    end
  #  end

  #  CMD.cmd("cut -f 1,#{good_field_pos*","} | uniq", :in => Misc.sort_stream(io), :pipe => true)
  #end

  #dep :mutation_basic_info
  #dep :mutation_incidence
  #task :mutation_info => :tsv do
  #  TSV.paste_streams(dependencies, :sort => true)
  #end

  #dep :mutation_info 
  #task :gene_status => :tsv do
  #  TSV.traverse mutation_info do |mutation, values|
  #    values = values.dup
  #    samples = values.pop

  #    new_values = []
  #    Misc.zip_fields(values)
  #    [gene, new_values]
  #  end
  #end


  #dep do |jobname,options|
  #  study = Study.setup(jobname.dup)
  #  jobs = study.genotyped_samples.collect{|sample| Sample.setup(sample, :cohort => study); sample.mutation_genes(:job) }.flatten
  #  Misc.bootstrap(jobs) do |job|
  #    job.run(true) unless (job.done? or job.started?) and not job.dirty?
  #    job.join unless job.done?
  #  end
  #  jobs
  #end
  #task :mutation_info_old => :tsv do

  #  parser = TSV::Parser.new dependencies.first
  #  fields = parser.fields - ["homozygous"]
  #  fields + ["Sample"]
  #  dumper = TSV::Dumper.new parser.options.merge(:fields => fields)
  #  dumper.init
  #  TSV.traverse dependencies, :into => dumper, :bar => true do |job|
  #    sample = job.sample.split(":").last
  #    values = job.path.tsv(:unnamed => true).collect{|k,v| [k,v[0..-2]+[[sample]]]}.reject{|k,v| v[0].empty?}
  #    values.extend MultipleResult
  #    values
  #  end

  #  io = TSV.traverse TSV.collapse_stream(dumper.stream), :type => :array, :into => :stream do |line|
  #    if line =~ /^#/
  #      line
  #    else
  #      genes, *rest = line.split("\t").collect{|v| v.split("|")}
  #      samples = rest.pop
  #      new_parts = []
  #      seen = Set.new
  #      genes.each_with_index do |gene,i|
  #        next if seen.include? gene
  #        seen << gene
  #        new_parts << [gene] + rest.collect{|v| v[i] } + [samples]
  #      end
  #      Misc.zip_fields(new_parts).collect{|p| p *  "|"} * "\t"
  #    end
  #  end
  #end

  #dep :mutation_info
  #task :gene_mutation_status => :tsv do
  #  mutation_info = step(:mutation_info)
  #  fields = TSV.parse_header(mutation_info).fields - ["Ensembl Gene ID"]
  #  dumper = TSV::Dumper.new :key_field => "Ensembl Gene ID", :fields => fields + ["Sample"], :type => :double, :namespace => organism
  #  dumper.init
  #  TSV.traverse mutation_info, :into => dumper, :bar => true do |mutation,values|
  #    values = values.dup
  #    samples = values.pop
  #    res = []
  #    Misc.zip_fields(values).each do |gene, *rest|
  #    end
  #    res
  #  end
  #  TSV.collapse_stream dumper.stream
  #end

  #dep do |jobname,options|
  #  study = Study.setup(jobname.dup)
  #  jobs = study.genotyped_samples.collect{|sample| Sample.setup(sample, :cohort => study); sample.mutation_genes(:job) }.flatten
  #  jobs
  #end
  #task :gene_mutation_status_old => :tsv do
  #  fields = dependencies.first.path.tsv.fields - ["Ensembl Gene ID"]
  #  dumper = TSV::Dumper.new :key_field => "Ensembl Gene ID", :fields => fields + ["Sample"], :type => :double, :namespace => organism
  #  dumper.init
  #  TSV.traverse dependencies, :into => dumper, :bar => true do |job|
  #    sample = job.sample
  #    res = []
  #    res.extend MultipleResult
  #    TSV.traverse job.path, :key_field => "Ensembl Gene ID", :fields => fields, :bar => true do |genes,lvalues|
  #      genes.each_with_index do |gene,i|
  #        res << [gene, lvalues.collect{|v| [v[i]]} + [[sample.split(":").last]]]
  #      end
  #    end
  #    res
  #  end
  #  TSV.collapse_stream dumper.stream
  #end

  #dep :gene_mutation_status 
  #task :sample_genes1 => :tsv do
  #  io = TSV.reorder_stream TSV.get_stream(step(:gene_mutation_status)), {6 => 0}
  #  parser = TSV::Parser.new io
  #  dumper = TSV::Dumper.new parser.options
  #  dumper.init
  #  TSV.traverse parser, :into => dumper do |ks,values|
  #    res = []
  #    res.extend MultipleResult
  #    ks.each_with_index do |k,i|
  #      res << [k, values.collect{|v| v.length != 1 ? v[i] : v.first }]
  #    end
  #    res
  #  end
  #  io2 = TSV.reorder_stream dumper.stream, {6 => 1}
  #  TSV.collapse_stream io2
  #end

  #dep :sample_genes1
  #dep do |jobname,options|
  #  study = Study.setup(jobname.dup)
  #  jobs = study.genotyped_samples.collect{|sample| Sample.setup(sample, :cohort => study); [sample.broken_genes(:job), sample.surely_broken_genes(:job)] }.flatten
  #  #Misc.bootstrap(jobs, 20, :bar => "Boostrapping jobs for sample_genes1") do |job|
  #  #  job.run(true) unless job.done? and not job.dirty? and not job.started?
  #  #end
  #  #Step.wait_for_jobs jobs
  #  jobs
  #end
  #task :sample_genes2 => :tsv do 
  #  Step.wait_for_jobs dependencies
  #  jobs = dependencies.dup
  #  mutation_status = dependencies.shift

  #  broken_genes = {}
  #  surely_broken_genes = {}
  #  dependencies.each do |dep|
  #    sample = dep.sample.split(":").last
  #    case dep.task.name
  #    when :broken_genes
  #      broken_genes[sample] = dep.load
  #    when :surely_broken_genes
  #      surely_broken_genes[sample] = dep.load
  #    else
  #      raise dep.task.name.to_s + ' not recognized'
  #    end
  #  end
  #  parser = TSV::Parser.new mutation_status
  #  fields = parser.fields + ["missing?", "missing"]
  #  dumper = TSV::Dumper.new parser.options.merge(:fields => fields)
  #  dumper.init
  #  TSV.traverse mutation_status, :into => dumper do |sample, values|
  #    sample = sample.first if Array === sample
  #    genes = values.first

  #    broken = broken_genes[sample]
  #    surely_broken = surely_broken_genes[sample]
  #    missing = genes.collect{|g| broken.include? g}
  #    smissing = genes.collect{|g| surely_broken.include? g}

  #    [sample, values + [missing, smissing]]
  #  end
  #end

  #dep :sample_genes2
  #dep :TSS
  #dep :mutation_incidence
  #task :sample_genes => :tsv do 
  #  Step.wait_for_jobs dependencies

  #  jsample_genes, jtss, jincidence = dependencies

  #  tss = jtss.load

  #  index = Association.index jincidence.path, {}, :file => jincidence.file(:index)

  #  parser = TSV::Parser.new jsample_genes
  #  fields = parser.fields + ["TSS"]
  #  dumper = TSV::Dumper.new parser.options.merge(:fields => fields)
  #  dumper.init
  #  TSV.traverse parser, :into => dumper do |sample, values|
  #    sample = sample.first if Array === sample
  #    genes = values.first
  #    mutations = index.reverse.match(sample).collect{|p| p.partition("~").last}
  #    tss_genes = Set.new tss.values_at(*mutations)
  #    tss_values = genes.collect{|g| tss_genes.include? g}
  #    values = values + [tss_values]

  #    missing = tss_genes - genes

  #    if missing.any?
  #      missing.each do |gene|
  #        values[0] << gene
  #        values[1..-2].each do |list|
  #          list << false
  #        end
  #        values[-1] << true
  #      end
  #    end
  #    [sample, values]
  #  end
  #end

end
