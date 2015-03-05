module Study

  input :database, :string, "Database code", "Appris"
  input :feature, :string, "Feature name", "feature"
  task :feature_report => :tsv do |database,feature|
    tsv = TSV.setup({}, :key_field => "Mutated Isoform", :fields => ["Position", "Description", "Sample"], :type => :double, :namespace => organism)

    Misc.bootstrap study.genotyped_samples do |sample|
      sample.structure_annotations(:principal => true)
      sample.structure_neighbour_annotations(:principal => true)
    end

    TSV.traverse study.genotyped_samples, :type => :array, :bar => "Extracting #{[feature, database]*"@"} annotations" do |sample|

      matches = {}

      sa = sample.structure_annotations(:principal => true)
      sa_fields = sa.fields
      organism = sa.namespace
      db_fields = sa.fields.select{|f| f =~ /#{database}/}.collect{|f| sa.fields.index f}
      sa.each do |mi,annotations|
        annotations = annotations.values_at *db_fields
        zipped = Misc.zip_fields(annotations)
        zipped.each do |type,pos,desc|
          next unless type == feature
          matches[mi] ||= []
          matches[mi] << [pos, desc]
        end
      end

      nmatches = {}
      sa = sample.structure_neighbour_annotations(:principal => true)
      sa_fields = sa.fields
      db_fields = sa.fields.select{|f| f =~ /#{database}/}.collect{|f| sa.fields.index f}
      sa.each do |mi,annotations|
        annotations = annotations[0..3]
        annotations = annotations.values_at *([0] + db_fields)
        zipped = Misc.zip_fields(annotations)
        zipped.each do |res,ltype,lpos,ldesc|
          Misc.zip_fields([ltype.split(";",-1), lpos.split(";",-1), ldesc.split(";",-1)]).each do |type,pos,desc|
            next unless type == feature
            nmatches[mi] ||= []
            nmatches[mi] << [res, pos, desc]
          end
        end
      end

      matches.each do |mi, vs|
        mi = MutatedIsoform.setup(mi.dup, :organism => organism)
        vs.each do |v|
          tsv.zip_new mi, v + [sample]
        end
      end

      nmatches.each do |mi, vs|
        mi = MutatedIsoform.setup(mi.dup, :organism => organism)
        vs.each do |v|
          tsv.zip_new mi, v[1..2] + [sample]
        end
      end
    end

    tsv
  end

  dep :feature_report
  task :feature_summary => :tsv do
    report = step(:feature_report).load
    summary = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => ["Feature", "Sample"], :namespace => organism, :type => :double)

    TSV.traverse report, :into => summary do |mi,values|
      pos, desc, sample = values
      gene = mi.protein.gene
      next if gene.nil?
      [gene, [desc,sample]]
    end

    summary.attach study.sample_info.to_double
  end

  dep :feature_summary
  extension :svg
  task :feature_plot => :text do
    summary = step(:feature_summary).load
    matrix = TSV.setup({}, :key_field => "Key", :fields => ["Ensembl Gene ID", "Feature", "Condition", "Sample"], :type => :list, :namepace => summary.namespace)

    summary = summary.select{|g,values| values["Sample"].uniq.length > 1}

    TSV.traverse summary, :into => matrix do |gene,lvalues|
      result = []
      result.extend MultipleResult
      seen = []
      Misc.zip_fields(lvalues).each do |values|
        feature, sample, condition = values
        next if seen.include? sample
        seen << sample
        key = [gene, feature, condition, sample] * ":"
        result  << [key, [gene, feature, condition, sample]]
      end
      result
    end


    matrix.attach Organism.identifiers(summary.namespace), :fields => ["Associated Gene Name"]

    matrix.R <<-EOF

library(ggplot2)

colnames(data) <- make.names(colnames(data))
plot <- ggplot(data=data, aes(x=Associated.Gene.Name, fill=Condition,color=Feature)) + geom_bar()
plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(plot=plot, filename='#{path}')

    EOF

    nil
  end

  task :kinmut => :tsv do
    Misc.bootstrap study.genotyped_samples do |sample|
      sample.ns_mutated_isoforms(:principal => true)
    end

    Workflow.require_workflow "KinMut"
    mis = []
    study.genotyped_samples.each do |sample|
      mis.concat  sample.ns_mutated_isoforms(:principal => true)
    end
    mis.uniq!

    predictions = KinMut.job(:predict_tsv, study, :list => mis).run

    result = TSV.setup({}, :key_field => "Mutated Isoform", :fields => ["Sample", "Damaged", "Prediction"], :type => :double, :namespace => organism)
    study.genotyped_samples.each do |sample|
      mis = sample.ns_mutated_isoforms(:principal => true)
      sample_preds = {}
      mis.each{|mi| sample_preds[mi] = predictions[mi] if predictions[mi]}
      sample_preds.each do |mi, pred|
        result.zip_new mi, [sample, pred.to_f >= 0 ? "true" : "false", pred]
      end
    end


    result
  end

  dep :kinmut
  task :kinmut_summary => :tsv do
    kinmut = step(:kinmut).load

    summary = TSV.setup({}, :key_field => "Ensembl Gene ID", :fields => ["Sample", "Mutated Isoform", "Damaged"], :namespace => organism, :type => :double)

    TSV.traverse kinmut, :into => summary do |mi,lvalues|
      res = []
      Misc.zip_fields(lvalues).each do |values|
        sample, damaged = values
        gene = mi.protein.gene
        next if gene.nil?
        res << [gene, [sample, mi, damaged]]
      end
      res.extend MultipleResult
      res
    end

    summary.attach study.sample_info.to_double
  end

  dep :kinmut_summary
  extension :svg
  input :damaged, :boolean, "Damaged only", false
  task :kinmut_plot => :text do |damaged|
    summary = step(:kinmut_summary).load
    matrix = TSV.setup({}, :key_field => "Key", :fields => ["Ensembl Gene ID", "Feature", "Condition", "Sample"], :type => :list, :namepace => summary.namespace)

    summary = summary.select{|g,values| values["Sample"].uniq.length > 10} if not damaged

    TSV.traverse summary, :into => matrix do |gene,lvalues|
      result = []
      result.extend MultipleResult
      seen = []
      Misc.zip_fields(lvalues).each do |values|
        sample, mi, feature, condition = values
        next if seen.include? sample
        next if damaged and feature == 'false'
        seen << sample
        key = [gene, feature, condition, sample] * ":"
        result  << [key, [gene, feature, condition, sample]]
      end
      result
    end


    matrix.attach Organism.identifiers(summary.namespace), :fields => ["Associated Gene Name"]

    matrix.R <<-EOF

library(ggplot2)

colnames(data) <- make.names(colnames(data))
plot <- ggplot(data=data, aes(x=Associated.Gene.Name, fill=Condition,opacity=Feature)) + geom_bar()
plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(plot=plot, filename='#{path}')

    EOF

    nil
  end

  task :broken_ppi => :tsv do 
    Misc.bootstrap study.genotyped_samples do |sample,i|
      job = sample.interfaces(:job, :principal => true)
      next if job.path.exists?
      job.run(false) 
    end

    tsv = TSV.setup({}, :key_field => "Sample", :fields => ["Ensembl Protein ID", "Partner Ensembl Protein ID"], :type => :double, :namespace => organism) 

    study.genotyped_samples.each do |sample|
      job = sample.interfaces(:job, :principal => true)
      TSV.traverse job.path, :fields => ["Partner Ensembl Protein ID"], :type => :flat, :into => tsv do |mi,partners|
        protein, _sep, change = mi.partition ":"
        next if protein =~ /ENST/
        pairs = Misc.zip_fields(partners.collect{|p| [protein, p]}.uniq)
        [sample, pairs]
      end
    end
    tsv
  end

  dep :broken_ppi
  extension :svg
  task :broken_ppi_plot => :binary do 
    summary = step(:broken_ppi).load
    matrix = TSV.setup({}, :key_field => "Key", :fields => ["Gene", "Partner", "Sample"], :type => :list, :namepace => summary.namespace)

    p2g = Organism.identifiers(organism).index :persist => true, :target => "Associated Gene Name", :fields => ["Ensembl Protein ID"]

    good = Misc.counts(summary.values.collect{|v| v.flatten.uniq}.flatten).select{|g,c| c > 6}.collect{|g,c| g}

    TSV.traverse summary, :into => matrix do |sample,lvalues|
      result = []
      result.extend MultipleResult
      Misc.zip_fields(lvalues).each do |values|
        p, pp = values
        next unless good.include? p
        gene = p2g[p]
        partner = p2g[pp]
        key = [gene, partner, sample] * ":"
        result  << [key, [gene, partner, sample]]
      end
      result
    end

    matrix = matrix.attach study.sample_info.to_double

    matrix.R <<-EOF

library(ggplot2)

colnames(data) <- make.names(colnames(data))
data = rbind(data, data.frame(Gene=data$Partner, Partner=data$Gene, Sample=data$Sample, Condition=data$Condition))
plot <- ggplot(data=data, aes(x=Gene, fill=Condition)) + geom_bar()
plot <- plot + theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(plot=plot, filename='#{path}', width=12, height=10)

    EOF

    nil
  end


end
