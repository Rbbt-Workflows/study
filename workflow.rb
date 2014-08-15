require 'rbbt-util'
require 'rbbt/workflow'

require 'study'
module Study
  extend Workflow

  helper :study do
    Study.setup(clean_name.dup)
  end

  helper :study_info do
    Study.study_info(study)
  end

  helper :organism do
    Study.organism(study)
  end

  helper :samples do
    Study.samples(study)
  end
end

require 'study/tasks/genotypes'
Workflow.require_workflow "MutationSignatures"
require 'study/tasks/mutation_signatures'

module Study
  extend Entity 

  property :job do |name|
    Study.job(name.to_sym, self)
  end

  Study.tasks.each do |name, b|
    property name.to_sym => :single do |run=true|
      job = job(name)
      case run
      when nil, TrueClass
        job.run 
      when :path
        job.run(true).join.path
      when :job
        job
      end
    end
  end
end

require 'rbbt/entity/study'


