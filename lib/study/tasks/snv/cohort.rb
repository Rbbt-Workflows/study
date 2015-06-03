CohortTasks = Proc.new do

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
  dep Enrichment, :rank_enrichment, :organism => :organism, :list => :sorted_significant_genes, :background => :mappable_genes
  input :database, :string, "Database to use", nil, :select_options => Enrichment::DATABASES
  task :significance_rank_enrichment => :tsv do
    TSV.get_stream step(:rank_enrichment)
  end

end
