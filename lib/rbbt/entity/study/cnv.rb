module Study
  def cnv_samples
    samples.select{|s| s.has_cnv?}
  end

  def has_cnv?
    cnv_samples.any?
  end
end
