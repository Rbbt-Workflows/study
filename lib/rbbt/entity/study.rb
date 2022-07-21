require 'rbbt/entity'
Workflow.require_workflow "Sample"
require 'rbbt/entity/sample'

$study_kb_dir ||= SOPT.get("-skbd--study_kb_dir* Directory where study knowledge-bases are stored ")[:study_kb_dir]
$study_kb_dir = Path.setup($study_kb_dir) if $study_kb_dir

module Study
  extend Entity

  attr_accessor :knowledge_base

  class << self 
    attr_accessor :knowledge_base_dir
  end
  self.knowledge_base_dir = $study_kb_dir || Rbbt.var.knowledge_base.studies

  def self.setup(obj, options = {})
    obj.extend Study
    obj
  end

  def matrices(*args)
    Study.matrices(self)
  end

  def has_matrix?(name)
    matrices.include? name
  end

  def matrix(matrix, *rest)
    Study.matrix(self, matrix, *rest)
  end

  def match_samples(list, target_list, through_field = nil)
    return list if (list & target_list).any?
    source_field, count = TSV.guess_id sample_info, list
    target_field, count = TSV.guess_id sample_info, target_list
    if through_field
      index_1 = sample_info.reorder(source_field, [through_field], :merge => true).to_flat
      index_2 = sample_info.reorder(through_field, [target_field], :merge => true).to_flat

      tmp = index_1.chunked_values_at(list).compact.flatten
      res = index_2.chunked_values_at(tmp).compact.flatten
      res
    else
      sample_info.index(:target => target_field).values_at *list
    end
  end

  def knowledge_base
    @@knowledge_bases ||= {}
    @knowledge_base ||= 
      @@knowledge_bases[self] ||= begin
                                    knowledge_base = KnowledgeBase.new Study.knowledge_base_dir[self], self.organism

                                    knowledge_base.format["Gene"] = "Ensembl Gene ID"
                                    knowledge_base.entity_options["Sample"] = {:cohort => self, :organism => self.organism }
                                    knowledge_base.entity_options["Gene"] = {:cohort => self, :organism => self.organism }
                                    knowledge_base.entity_options["GenomicMutation"] = {:watson => watson, :organism => self.organism }

                                    if self.has_genotypes?
                                      
                                      knowledge_base.register :mutation_info, nil, :source => "Genomic Mutation", :taget => "Ensembl Gene ID", :merge => false do 
                                        d = TSV.reorder_stream_tsv self.job(:mutation_info).produce, "Genomic Mutation"
                                        TSV.collapse_stream(d)
                                      end

                                      knowledge_base.register :sample_mutations, nil, :source => "Sample", :taget => "Genomic Mutation", :merge => false do 
                                        d = TSV.reorder_stream_tsv self.job(:mutation_incidence).produce, "Sample"
                                        TSV.collapse_stream(d)
                                      end

                                      knowledge_base.register :sample_genes, nil, :source => "Sample", :target => "Gene", :merge => false do 
                                        d = TSV.reorder_stream_tsv self.job(:sample_genes).produce, "Sample"
                                        TSV.collapse_stream(d)
                                      end
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

  property :sample_extended_info => :single do
    Study.sample_extended_info(self)
  end


  property :metadata => :single do
    study_info
  end

  property :users => :single do
    Study.users(self)
  end

  property :samples => :single do
    samples = Sample.setup(Study.samples(self), self)
    samples.extend AnnotatedArray
    samples
  end

  property :organism => :single do
    Study.organism(self)
  end

  property :subset => :single do |database,options={}|
    options = Misc.add_defaults options, :source => :all, :target => :all
    self.knowledge_base.subset(database, options)
  end

  property :mutations_in_range => :single do |chr,start,eend|
    TSV.traverse self.genomic_mutations, :into => [], :type => :array do |mutation|
      mut_chr,pos,allele = mutation.split(":")
      next unless mut_chr == chr
      pos = pos.to_i
      next unless pos >= start and pos <= eend
      mutation 
    end
  end

end

require 'rbbt/entity/study/genotypes'
require 'rbbt/entity/study/cnv'

if defined? Sample and Entity === Sample
  module Sample
    self.format = ["submitted_specimen_id", "submitted_sample_id", "icgc_specimen_id", "icgc_donor_id"]

    property :study => :both do
      Study.setup(cohort)
    end
  end
end


if defined? Gene and Entity === Gene
  module Gene

    property :significant_in_study => :array do |study|
      sg = Study.setup(study).significant_genes
      self.collect do |gene|
        sg.include? gene
      end
    end

    property :damage_bias_in_study => :array do |study|
      db = Study.setup(study).gene_damage_bias
      db.chunked_values_at(self).collect{|values| values.nil? ? nil : values.last}
    end

    property :in_study_property => :array2single do |study,property,*args|
      Study.setup(study.dup).send(property, *args).chunked_values_at self
    end
  end
end
