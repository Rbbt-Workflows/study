require 'rbbt/entity'
Workflow.require_workflow "Sample"
require 'rbbt/entity/sample'

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

  def matrices(*args)
    Study.matrices(self)
  end

  def matrix(*args)
    Study.matrix(self, *args)
  end

  def match_samples(list)
    sample_info.index(:target => "Sample").values_at *list
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
                                      #knowledge_base.register :mutation_info do
                                      #  job = self.job(:mutation_info)
                                      #  job.run(true) unless job.done? or job.started?
                                      #  Step.wait_for_jobs job unless job.done?
                                      #  fields = TSV.parse_header(job.path).fields - ["Ensembl Gene ID"]
                                      #  TSV.open job.path, :key_field => "Genomic Mutation", :fields => ["Ensembl Gene ID"] + fields
                                      #end

                                      knowledge_base.register :mutation_info, nil, :source => "Genomic Mutation", :taget => "Ensembl Gene ID" do self.job(:mutation_info) end

                                      #knowledge_base.register :sample_mutations_old do
                                      #  job = self.job(:mutation_incidence)
                                      #  job.run(true) unless job.done? or job.started?
                                      #  Step.wait_for_jobs job unless job.done?
                                      #  job.path.tsv :key_field => "Sample", :type => :flat
                                      #end

                                      knowledge_base.register :sample_mutations, nil, :source => "Sample", :taget => "Genomic Mutation" do self.job(:mutation_incidence) end

                                      #knowledge_base.register :sample_genes_old do
                                      #  job = self.job(:sample_genes)
                                      #  job.run(true) unless job.done? or job.started?
                                      #  Step.wait_for_jobs job unless job.done?
                                      #  job.path.tsv :key_field => "Sample", :merge => true
                                      #end
                                      knowledge_base.register :sample_genes, nil, :source => "Sample", :target => "Gene" do self.job(:sample_genes) end
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
