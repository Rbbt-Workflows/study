module Study
  dep do |jobname, task|
    study = Study.setup(jobname.dup)
    study.samples.collect{|sample| sample.gene_cnv_status(:job)}
  end
  task :sample_gene_cnv => :tsv do
    Step.wait_for_jobs dependencies
    gene_sample_cnv = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => ["Status", "Sample"], :type => :double, :namespace => organism)
                                
    TSV.traverse dependencies.collect{|dep| dep.path}, :into => gene_sample_cnv, :bar => self.progress_bar("Processing samples"), :cpus => 10 do |path|
      res = []
      sample = File.basename(path).sub(/_[^_]+$/,'').sub(/.*:/,'')
      TSV.traverse path, :into => res do |gene,statuses|
        next if statuses.nil? or statuses.empty?
        [gene, [statuses*";", sample]]
      end
      res.extend MultipleResult
      res
    end
  end

end
