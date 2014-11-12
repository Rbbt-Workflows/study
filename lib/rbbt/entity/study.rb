require 'rbbt/entity'

module Study
  extend Entity

  attr_accessor :knowledge_base

  class << self 
    attr_accessor :knowledge_base_dir
  end
  self.knowledge_base_dir = Rbbt.var.knowledge_base.studies

  def self.setup(obj, options = {})
    obj.extend Study
    obj
  end

  def knowledge_base
    @@knowledge_bases ||= {}
    @knowledge_base ||= 
      @@knowledge_bases[self] ||= begin
                                    knowledge_base = KnowledgeBase.new Study.knowledge_base_dir[self], self.organism

                                    knowledge_base.format["Gene"] = "Ensembl Gene ID"
                                    knowledge_base.entity_options["Sample"] = {:cohort => self }
                                    knowledge_base.entity_options["GenomicMutation"] = {:watson => watson }

                                    if self.has_genotypes?
                                      job = self.job(:mutation_info)
                                      job.run(true)

                                      fields = TSV.parse_header(job.path).fields - ["Sample"]
                                      knowledge_base.register :mutation_info, job.path, :source => "Genomic Mutation", :target => "Ensembl Gene ID", :fields => fields, :merge => true

                                      knowledge_base.register :sample_mutations, job.path, :source => "Sample", :target => "Genomic Mutation", :merge => true

                                      job = self.job(:sample_genes)
                                      job.run(true)
                                      
                                      knowledge_base.register :sample_genes, job.path
                                    end


                                    knowledge_base
                                  end
  end

  attr_accessor :dir
  def dir
    @dir ||= Study.study_dir[self].find
  end

  def job(name, options={})
    Study.job(name, self, options)
  end

  property :study_info => :single do
    Study.study_info(self)
  end

  property :sample_info => :single do
    Study.sample_info(self)
  end

  property :metadata => :single do
    study_info
  end

  property :users => :single do
    Study.users(self)
  end

  property :samples => :single do
    Sample.setup(Study.samples(self), self)
  end

  property :organism => :single do
    Study.organism(self)
  end

  property :subset => :single do |database,options={}|
    options = Misc.add_defaults options, :source => :all, :target => :all
  self.knowledge_base.subset(database, options)
  end

end

require 'rbbt/entity/study/genotypes'

if defined? Sample and Entity === Sample
  module Sample
    self.format = ["submitted_specimen_id", "submitted_sample_id", "icgc_specimen_id", "icgc_donor_id"]

    property :study => :both do
      Study.setup(cohort)
    end
  end
end
