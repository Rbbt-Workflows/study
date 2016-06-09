require 'sample/tasks/snv/common'
require 'study/tasks/snv/cohort'

module Study

  task :exome => :boolean do
    study = Study.setup(clean_name.dup)
    exome = study.metadata[:exome]
    exome = true if exome.nil? or (String === exome and exome.empty?)
    exome
    false
  end

  task :watson => :boolean do
    true
  end

  dep do |jobname, task|
    study = Study.setup(jobname.dup)
    study.genotyped_samples.collect{|sample| Sample.setup(sample, :cohort => study).genomic_mutations(:job)}
  end
  task :mutation_incidence => :tsv do
    
    size = Math.sqrt(dependencies.length).to_i + 1
    chunks = Misc.ordered_divide(dependencies, size)
    bar = self.progress_bar("Processing sample chunks", :max => chunks.length)
    TSV.traverse (0..chunks.length-1).to_a, :cpus => 5 do |i|
      bar.tick
      dep_list = chunks[i]
      streams = dep_list.collect do |dep|
        sample = File.basename(dep.path).split(":").last
        TSV.traverse CMD.cmd('env LC_ALL=C sort', :in => dep.path.open, :pipe => true), :into => :stream, :type => :array, :bar => "Processing #{ sample }" do |mutation|
          [mutation, sample] * "\t"
        end
      end

      io = TSV.paste_streams(streams, :sort => false, :key_field => "Genomic Mutation", :fields => ["Sample"], :same_fields => true, :type => :double, :merge => true, :namespace => organism)
      Open.write(file('tmp-stream-' + i.to_s), io.read)
    end

    log :pasting, "Pasting sorted streams"
    _tmp_streams = file('tmp-stream-*').glob.each{|f| f.open }
    TSV.paste_streams(_tmp_streams, :sort => false, :same_fields => true, :type => :double, :filename => path)
  end

  dep :mutation_incidence
  returns "Genomic Mutation"
  task :genomic_mutations => :array do 
    TSV.traverse step(:mutation_incidence), :type => :flat, :into => :stream do |mutation, samples|
      res = [mutation] * samples.length
      res.extend MultipleResult
      res
    end
  end

end

Study.instance_eval &SNVTasks
Study.instance_eval &CohortTasks
