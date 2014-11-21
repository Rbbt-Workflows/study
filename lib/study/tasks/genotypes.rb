module Study

  task :organism => :string do
    Study.organism(study)
  end

  task :metagenotype => :array do
    study.genotyped_samples.collect{|sample|
      begin
        sample = [study, sample] * ":" unless sample.include? study
        Sample.mutations(sample).read.split("\n")
      rescue Exception
        []
      end
    }.flatten
  end


  task :mutation_info => :tsv do
    tsv = nil
    
    sample_mutations = {}
    genotyped_samples = study.genotyped_samples
    raise "No genotyped samples" if genotyped_samples.empty?

    Misc.bootstrap genotyped_samples do |sample|
      sample.mutation_genes
    end

    genotyped_samples.each do |sample|
      sample_mutation_info = sample.mutation_genes
      sample_mutations[sample] = Set.new sample_mutation_info.keys
      good_fields = sample_mutation_info.fields - ["missing"]
      if tsv.nil?
        tsv = sample_mutation_info.slice(good_fields)
      else
        tsv.merge! sample_mutation_info.slice(good_fields)
      end
    end

    raise "No mutations" if tsv.nil?

    tsv.add_field "Sample" do |mutation, info|
      sample_mutations.select{|sample,mutations| mutations.include? mutation }.collect{|sample, mutations| sample }
    end 

    tsv
  end

  dep :mutation_info
  task :mutation_samples => :tsv do
    step(:mutation_info).join.path.tsv :fields => ["Sample"], :type => :flat
  end

  dep :mutation_info
  task :sample_genes => :tsv do
    # ToDo: Check when Sample is in list of fields
    
    log :loading, "Loading mutation info for #{ study }"

    fields = TSV.parse_header(step(:mutation_info).join.path).fields - ["Sample"]
    mutation_info = step(:mutation_info).join.path.tsv(:fields => fields, :unnamed => true)
    sample_mutations = step(:mutation_info).join.path.tsv(:key_field => "Sample", :fields => ["Genomic Mutation"], :merge => true, :zipped => true, :unnamed => true, :type => :double).to_flat

    fields << "missing?"
    fields << "missing"
    tsv = TSV.setup({}, :key_field => "Sample", :fields => fields, :namespace => organism, :type => :double)

    genotyped_samples = study.genotyped_samples

    TSV.traverse genotyped_samples, :into => tsv, :cpus => [genotyped_samples.length / 10, 5].min, :respawn => 1.8, :bar => "Sample genes for study" do |sample|
      Sample.setup sample, study

      mutations = broken = surely_broken = mutation_values_list = nil

      log :prepare, "Prepare info for #{ sample }" do
        mutations = sample_mutations[sample]
        broken = sample.broken_genes
        surely_broken = sample.surely_broken_genes
        mutation_values_list = mutation_info.select(:key => mutations).values
      end

      gene_mutations = {}
      values = nil
      
      log :gene_mutations, "Compute gene_mutation info for #{ sample }" do
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

    tsv
  end

  dep :organism
  dep :metagenotype
  dep Sequence, :binomial_significance, :organism => :organism, :mutations => :metagenotype 
  task :binomial_significance => :tsv do
    TSV.get_stream step(:binomial_significance)
  end

  Workflow.require_workflow "MutationEnrichment"
  dep :organism
  dep :mutation_samples
  dep MutationEnrichment, :sample_pathway_enrichment, :organism => :organism, :mutations => :mutation_samples 
  task :sample_enrichment => :tsv do
    TSV.get_stream step(:sample_pathway_enrichment)
  end
end
