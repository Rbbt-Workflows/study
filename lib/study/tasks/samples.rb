module Study

  task :organism => :string do
    Study.organism(study)
  end

  task :samples => :array do
    Study.samples(clean_name)
  end

  task :genotyped_samples => :array do
    Study.genotyped_samples(clean_name)
  end

  dep :genotyped_samples
  task :num_genotyped_samples => :integer do
    step(:genotyped_samples).load.length
  end
end
