require 'rbbt/matrix'
require 'rbbt/matrix/differential'

Workflow.require_workflow "Sample"

module Study

  class << self
    attr_accessor :study_dir

    def study_dir
      @study_dir ||= Sample.study_repo
    end

    def study_dir=(dir)
      dir = Path.setup dir.dup unless Path === dir
      Sample.study_repo = dir
      @study_dir = dir
    end
  end

  def self.find_study(study)
    return Path.setup(study) if File.exists?(study)
    return study_dir[study] if study_dir[study].exists?
    return Sample.project_repo["*"][study].glob.first if Sample.project_repo["*"][study].glob.any?
    raise "Study not found: #{ study }"
  end

  def self.study_info(study)
    path = find_study(study)
    path["metadata.yaml"].yaml
  end

  def self.users(study)
    study_info(study)[:users]
  end

  def self.genotyped_samples(study)
    path = find_study(study)
    samples = path.genotypes.vcf.glob("*").collect{|f| File.basename(f).sub(/\.vcf.*/,'') }
    return samples if samples.any?
    path.genotypes.glob("*").collect{|f| File.basename f }
  end

  def self.cnv_samples(study)
    path = find_study(study)
    path.CNV.glob("*").collect{|f| File.basename(f) }
  end

  def self.sample_info_file(study)
    find_study(study).samples
  end

  def self.sample_info(study)
    path = sample_info_file(study)
    organism = study_info(study)[:organism]

    if path.exists?
      tsv = path.tsv :namespace => organism 
      tsv.entity_options = {:cohort => study, :organism => organism}
      tsv
    else
      samples = genotyped_samples(study)
      tsv = TSV.setup(samples, :key_field => "Sample", :fields => [], :type => :list, :namespace => organism)
      tsv.entity_options = {:cohort => study, :organism => organism}
      tsv
    end
  end

  def self.sample_extended_info_files(study)
    find_study(study).glob('sample_*') + 
    find_study(study).glob('donor_*')  
  end

  def self.sample_extended_info(study)
    sample_extended_info = sample_info(study)
    return nil if sample_extended_info_files(study).nil?
    sample_extended_info_files(study).each do |file|
      sample_extended_info = sample_extended_info.attach file
    end
    sample_extended_info
  end

  def self.organism(study)
    study_info(study)[:organism]
  end

  def self.samples(study)
    sample_info(study).keys
  end

  def self.matrices(study)
    path = find_study(study)
    path.matrices.glob("*").collect{|f| File.basename(f)}
  end

  def self.matrix_file(study, matrix)
    path = find_study(study)
    path.matrices[matrix.to_s].data
  end

  def self.matrix(study, matrix, format = nil)
    file = matrix_file(study, matrix)
    sample_info = sample_info_file(study)
    value_type = study_info(study)[:expression_type]
    Matrix.new file.find, sample_info, value_type, format
  end

  def self.studies
    case study_dir
    when nil
      []
    when Path 
      study_dir.find_all.collect do |study_path| 
        study_path.glob('*').select{|f| f.directory? }
      end.flatten.collect{|f| study_dir.annotate f} 
    else
      Dir.glob(File.join(study_dir, "*"))
    end.sort.collect do |dir| 
      study = Study.setup(File.basename(dir))
      study.dir = study_dir.annotate(dir)
      study
    end
  end

end
