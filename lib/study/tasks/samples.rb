module Study

  task :organism => :string do
    Study.organism(study)
  end

  task :condition => :string do
    study.metadata[:condition]
  end

  task :num_genotyped_samples => :integer do
    study.genotyped_samples.length
  end
end
