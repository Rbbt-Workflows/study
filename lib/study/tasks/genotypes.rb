module Study

  task :mutation_info => :tsv do
    tsv = nil
    
    sample_mutations = {}
    genotyped_samples = study.samples.select{|s| s.has_genotype? }
    raise "No genotyped samples" if genotyped_samples.empty?

    Misc.bootstrap genotyped_samples do |sample|
      sample.mutation_genes
    end

    genotyped_samples.each do |sample|
      sample_mutation_info = sample.mutation_genes
      sample_mutations[sample] = Set.new sample_mutation_info.keys
      if tsv.nil?
        tsv = sample_mutation_info
      else
        tsv.merge! sample_mutation_info
      end
    end

    raise "No mutations" if tsv.nil?

    tsv.add_field "Sample" do |mutation, info|
      sample_mutations.select{|sample,mutations| mutations.include? mutation }.collect{|sample, mutations| sample }
    end 

    tsv
  end

  dep :mutation_info
  task :sample_genes => :tsv do
    # ToDo: Check when Sample is in list of fields
    fields = TSV.parse_header(step(:mutation_info).join.path).fields - ["Sample"]
    mutation_info = step(:mutation_info).join.path.tsv(:fields => fields)
    sample_mutations = step(:mutation_info).join.path.tsv(:key_field => "Sample", :fields => ["Genomic Mutation"], :merge => true, :unnamed => true, :type => :double).to_flat

    fields << "missing?"
    fields << "missing"
    tsv = TSV.setup({}, :key_field => "Sample", :fields => fields, :namespace => organism, :type => :double)

    genotyped_samples = study.samples.select{|s| s.has_genotype? }

    TSV.traverse genotyped_samples, :into => tsv, :cpus => 10 do |sample|
      Sample.setup sample
      mutations = sample_mutations[sample]
      broken = sample.broken_genes
      surely_broken = sample.surely_broken_genes
      mutation_values_list = mutation_info.select(:key => mutations).values
      gene_mutations = {}
      mutation_values_list.each do |mutation_values|
        Misc.zip_fields(mutation_values).each do |_v|
          gene, *values = _v
          gene_mutations[gene] ||= []
          gene_mutations[gene] << (values + [broken.include?(gene), surely_broken.include?(gene)])
        end
      end

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

      new_values = Misc.zip_fields(values)
      [sample,  new_values]
    end

    tsv
  end

  task :binomial_significance => :tsv do
    Sequence.job(:binomial_significance, study, :mutations => study.metagenotype.sort, :threshold => 0.1, :organism => organism).run
  end
end
