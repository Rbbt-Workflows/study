require 'study/tasks/pancancer'
module Study

  dep do |jobname, task|
    study = Study.setup(jobname.dup)
    study.genotyped_samples.collect{|sample| sample.genomic_mutations(:job)}
  end
  task :metagenotype => :array do
    TSV.traverse dependencies, :into => :stream do |dep|
      dep.path.read.strip
    end
  end


  dep do |jobname, task|
    study = Study.setup(jobname.dup)
    Misc.bootstrap(study.genotyped_samples, 20, :bar => "Boostrapping genotyped samples") do |sample|
      job = Sample.setup(sample, :cohort => jobname).mutation_genes(:job)
      job.run(true) unless job.done? and not job.dirty? and not job.started?
      job.join unless job.done? 
    end
    study.genotyped_samples.collect{|sample| sample.mutation_genes(:job)}
  end
  task :mutation_info => :tsv do
    tsv = nil
    
    genotyped_samples = study.genotyped_samples
    raise "No genotyped samples" if genotyped_samples.empty?

    tsv = genotyped_samples.first.mutation_genes.annotate({})
    tsv.type = :double
    tsv.unnamed = true
    good_fields = tsv.fields - ["missing"]

    sample_mutations = {}
    io = TSV.traverse genotyped_samples, :type => :array, :bar => "Mutation info", :into => :stream do |sample|
      sample_mutation_info_job = sample.mutation_genes(:job)
      sample_mutation_info_job.run(true).join
      sample_mutation_info = sample_mutation_info_job.path.tsv :fields => good_fields, :unnamed => true, :type => :double
      mutations = sample.genomic_mutations
      sample_mutations[sample] = Set.new Annotated.purge(mutations)

      res = []
      sample_mutation_info.each do |mutation,lvalues|
        res << Misc.zip_fields(lvalues).collect{|values| ([mutation] + values) * "\t"}
      end

      next if res.empty?
      res * "\n"
    end

    sorted = CMD.cmd('env LC_ALL=C sort -u', :in => io, :pipe => true)

    raise "No mutations" if tsv.nil?

    TSV.traverse sorted, :type => :array, :bar => "Adding to TSV", :into => tsv do |line|
      next if line.empty?
      mut, *values = line.split("\t")
      _v = values.collect{|v| v.split("|")}
      [mut, _v]
    end

    tsv.with_monitor do
      tsv.add_field "Sample" do |mutation, info|
        sample_mutations.select{|sample,mutations| mutations.include? mutation }.collect{|sample, mutations| sample }
      end 
    end

    tsv
  end

  dep :mutation_info
  task :mutation_samples => :tsv do
    step(:mutation_info).join.path.tsv :fields => ["Sample"], :type => :flat, :unnamed => true
  end

  dep :mutation_info
  task :sample_genes => :tsv do
    mutation_info_job = step(:mutation_info)
    Step.wait_for_jobs mutation_info_job

    # ToDo: Check when Sample is in list of fields
    fields = TSV.parse_header(mutation_info_job.path).fields - ["Sample"]

    log :loading, "Loading sample mutations for #{ study }"
    sample_mutations = mutation_info_job.path.tsv(:key_field => "Sample", :fields => ["Genomic Mutation"], :merge => true, :zipped => true, :unnamed => true, :type => :double, :monitor => true).to_flat

    log :loading, "Loading mutation info for #{ study }"
    mutation_info = mutation_info_job.path.tsv(:fields => fields, :unnamed => true, :monitor => true)

    fields << "missing?"
    fields << "missing"

    genotyped_samples = study.genotyped_samples

    TSV.traverse genotyped_samples, :into => :dumper, 
      :key_field => "Sample", :fields => fields, :namespace => organism, :type => :double,
      :bar => self.progress_bar("Sample genes for #{ study }") do |sample|

      log [study, sample] * ": "
      Sample.setup sample, study

      mutations = broken = surely_broken = mutation_values_list = nil

      log :prepare, "Prepare info for #{ sample }" do
        mutations = sample_mutations[sample]
        broken = Annotated.purge(sample.broken_genes)
        surely_broken = Annotated.purge(sample.surely_broken_genes)
        mutation_values_list = mutation_info.select(:key => Annotated.purge(mutations)).values
      end

      gene_mutations = {}
      values = nil
      
      log :gene_mutations, "Compute gene_mutation info for #{ sample }: #{mutation_values_list.length}" do
        mutation_values_list.each do |mutation_values|
          Misc.zip_fields(mutation_values).each do |_v|
            gene, *values = _v
            gene_mutations[gene] ||= []
            gene_mutations[gene] << (values + [broken.include?(gene), surely_broken.include?(gene)])
          end
        end
      end

      log :gene_values, "Compute gene based values for #{ sample }" do
        values = gene_mutations.collect{|gene,list| 
          first = list.first.dup
          merged = list.inject(first) do |acc,e|
            e.each_with_index do |c,i|
              acc[i] = true if (acc[i] and acc[i] != "false") or c == "true"
            end
            acc
          end
          merged.unshift gene
          merged.flatten
        }
      end

      new_values = Misc.zip_fields(values)

      [sample,  new_values]
    end
  end

  #dep :organism
  #dep :metagenotype
  #dep :exome
  #dep Sequence, :binomial_significance, :organism => :organism, :mutations => :metagenotype, :exome => :exome
  #task :binomial_significance => :tsv do
  #  TSV.get_stream step(:binomial_significance)
  #end

  #Workflow.require_workflow "MutationEnrichment"
  #dep :organism
  #dep :mutation_incidence
  #dep MutationEnrichment, :sample_pathway_enrichment, :organism => :organism, :mutations => :mutation_incidence
  #input :database, :string, "Database to use", nil, :select_options => MutationEnrichment::DATABASES
  #task :sample_enrichment => :tsv do
  #  TSV.get_stream step(:sample_pathway_enrichment)
  #end
end
