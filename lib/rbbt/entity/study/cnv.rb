module Study
  def has_cnv?
    samples.select{|s| s.has_cnv?}.any?
  end
end
