
module Study
  def has_expression?
    Study.matrix_file(self, :gene_expression)
  end
end
