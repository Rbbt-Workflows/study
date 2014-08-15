Workflow.require_workflow "Sample"

module Study

  class << self
    attr_accessor :study_dir

    def study_dir
      @study_dir ||= Sample.study_repo
    end

    def study_dir=(dir)
      dir = Path.setup dir.dup unless Path === dir
      @study_dir = dir
    end
  end

  def self.find_study(study)
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

  def self.sample_info(study)
    path = find_study(study)
    organism = study_info(study)[:organism]

    return path.samples.tsv :namespace => organism if path.samples.exists?
    samples = path.genotypes.glob("*").collect{|f| File.basename f }
    tsv = TSV.setup(samples, :key_field => "Sample", :fields => [], :type => :list, :namespace => organism)
    tsv.entity_options = {:cohort => study, :organism => organism}
    tsv
  end

  def self.organism(study)
    study_info(study)[:organism]
  end

  def self.samples(study)
    sample_info(study).keys
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
