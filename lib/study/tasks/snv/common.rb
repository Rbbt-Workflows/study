Workflow.require_workflow "Sequence"
Workflow.require_workflow "Structure"
Workflow.require_workflow "GERP"
Workflow.require_workflow "DbSNP"
Workflow.require_workflow "DbNSFP"
Workflow.require_workflow "EVS"
SNVTasks = Proc.new do

  dep :genomic_mutations
  dep :organism 
  dep GERP, :annotate, :mutations => :genomic_mutations, :organism => :organism
  task :annotate_GERP => :tsv do
    TSV.get_stream step(:annotate)
  end
  
  dep :genomic_mutations
  dep :organism 
  dep DbSNP, :annotate, :mutations => :genomic_mutations, :organism => :organism
  task :annotate_DbSNP => :tsv do
    TSV.get_stream step(:annotate)
  end

  dep :genomic_mutations
  dep :organism 
  dep Genomes1000, :annotate, :mutations => :genomic_mutations, :organism => :organism
  task :annotate_Genomes1000 => :tsv do
    TSV.get_stream step(:annotate)
  end

  dep :genomic_mutations
  dep :organism 
  dep EVS, :annotate, :mutations => :genomic_mutations, :organism => :organism
  task :annotate_EVS => :tsv do
    TSV.get_stream step(:annotate)
  end

  dep :annotate_DbSNP
  dep :annotate_Genomes1000
  dep :annotate_GERP
  dep :annotate_EVS
  task :genomic_mutation_annotations => :tsv do
    TSV.paste_streams dependencies, :sort => true
  end

  dep :genomic_mutations
  dep Sequence, :mutated_isoforms_fast, :mutations => :genomic_mutations, :organism => :organism
  task :genomic_mutation_consequence => :tsv do
    TSV.get_stream step(:mutated_isoforms_fast)
  end

  dep :genomic_mutation_consequence
  task :mi => :array do
    io = TSV.traverse step(:genomic_mutation_consequence), :into => :stream do |mut, mis|
      mis = mis.reject{|mi| mi =~ /:([A-Z])\d+\1$/}
      mis.extend MultipleResult
      mis
    end
    CMD.cmd('shuf', :in => io, :pipe => true)
  end

  dep :mi
  dep :organism
  dep DbNSFP, :annotate, :mutations => :mi, :organism => :organism
  task :DbNSFP => :tsv do
    TSV.get_stream step(:annotate)
  end

  Workflow.require_workflow "KinMut_2"
  dep :mi
  dep KinMut2, :predict_fix, :mutations => :mi
  task :kinmut => :tsv do
    TSV.get_stream step(:predict_fix)
  end

  dep :DbNSFP
  input :field, :string, "Damage score field from DbNSFP", "MetaSVM_score"
  task :mi_damaged => :array do |field|
    TSV.traverse step(:DbNSFP), :fields => [field], :type => :single, :cast => :to_f, :into => :stream do |mi, score|
      next nil unless score > 0
      mi.extend MultipleResult if Array === mi
      mi
    end
  end

  dep :mi
  dep :organism
  dep Structure, :mi_interfaces, :mutated_isoforms => :mi, :organism => :organism
  task :interfaces => :tsv do
    parser = TSV::Parser.new step(:mi_interfaces)
    dumper = TSV::Dumper.new parser.options.merge(:fields => ["Partner Ensembl Protein ID"])
    dumper.init
    TSV.traverse parser, :into => dumper do |mi, values|
      mi = mi.first if Array === mi
      [mi, [values[1].uniq]]
    end
  end

  dep :mi
  dep :organism
  dep Structure, :annotate_mi, :mutated_isoforms => :mi, :organism => :organism, :database => "Appris"
  task :firestar => :tsv do
    fields = ["Appris Feature", "Appris Feature Description", "Appris Feature Range"]
    parser = TSV::Parser.new step(:annotate_mi)
    dumper = TSV::Dumper.new parser.options
    dumper = TSV::Dumper.new parser.options.merge(:fields => ["Firestar site", "Firestar range"])
    dumper.init
    TSV.traverse parser, :fields => fields, :into => dumper do |mi, values|
      next unless values[0].include? "firestar"
      mi = mi.first if Array === mi

      filtered = []
      Misc.zip_fields(values).each do |name,range,desc|
        next unless name == 'firestar'
        filtered << [desc, range]
      end

      next if filtered.empty?
      [mi, Misc.zip_fields(filtered)]
    end
  end

  dep :mi
  dep :organism
  dep Structure, :annotate_mi_neighbours, :mutated_isoforms => :mi, :organism => :organism, :database => "Appris"
  task :firestar_neighbours => :tsv do
    fields = ["Appris Feature", "Appris Feature Description", "Appris Feature Range"]
    parser = TSV::Parser.new step(:annotate_mi_neighbours)
    dumper = TSV::Dumper.new parser.options.merge(:fields => ["Firestar neighbour site", "Firestar neighbour range"])
    dumper.init
    TSV.traverse parser, :fields => fields, :into => dumper do |mi, values|
      next unless values[1].include? "firestar"
      mi = mi.first if Array === mi

      filtered = []
      Misc.zip_fields(values).each do |res,name,range,desc|
        next unless name == 'firestar'
        filtered << [desc, range]
      end

      next if filtered.empty?
      [mi, Misc.zip_fields(filtered)]
    end
  end

  dep :interfaces
  dep :firestar
  dep :firestar_neighbours
  task :mi_annotations => :tsv do
    TSV.paste_streams dependencies, :sort => true
  end

  dep :genomic_mutations
  dep :organism
  dep Sequence, :TSS, :positions => :genomic_mutations, :organism => :organism
  task :TSS => :tsv do
    TSV.get_stream step(:TSS)
  end
end