require 'rbbt-util'
require 'rbbt/workflow'


Workflow.require_workflow "Sample"

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

require 'study/tasks/samples'
#require 'study/tasks/genotypes'
#require 'study/tasks/cnv'
#Workflow.require_workflow "MutationSignatures"
#require 'study/tasks/mutation_signatures'
#require 'study/tasks/expression'

require 'study/tasks/snv'
require 'study/tasks/maf'

module Study
  extend Entity 

  property :job do |name,options={}|
    Study.job(name.to_sym, self,options)
  end

  def self.update_task_properties
    Study.tasks.each do |name, b|
      property name.to_sym => :single do |run=true, options={}|
        run, options = true, run if Hash === run

        job = job(name,options)
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
end

Study.update_task_properties

require 'rbbt/entity/study'


