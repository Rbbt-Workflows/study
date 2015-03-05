require 'rbbt/matrix'
require 'rbbt/matrix/barcode'
module Study

  task :barcode => :tsv do
    return nil unless Study.study_dir.has_expression?
    Study.matrix(study, :gene_expression).barcode(path)
    nil
  end
end
