$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), 'lib'))
$LOAD_PATH.unshift(File.dirname(__FILE__))
require 'test/unit'
require 'rbbt/workflow'

Workflow.require_workflow "Study"
require 'rbbt/entity/study'

class TestStudy < Test::Unit::TestCase
  def _test_bladder_sample_genes
    study = Study.setup("bladder-preal-17samples")
    assert study.knowledge_base.get_index(:sample_genes, :persist => false).include? "11401~ENSG00000181555"

    index = study.knowledge_base.get_index(:sample_genes, :persist => false)
    index.unnamed = false
    assert_equal "false", index["11401~ENSG00000181555"]["affected"]
    assert_equal "true", index["11401~ENSG00000183337"]["affected"]
  end

  def _test_bladder_sample_genes_names
    study = Study.setup("bladder-preal-17samples")
    kb = study.knowledge_base
    assert kb.get_index(:sample_genes, :persist => false, :target_format => "Associated Gene Name").include? "11401~BCOR"
  end

  def test_CLL
    study = Study.setup("Chronic_Lymphocytic_Leukemia-ISC_MICINN-ES")
    kb = study.knowledge_base
    index = kb.get_index(:sample_genes, :persist => false, :target_format => "Associated Gene Name")
    Log.tsv index
    sf3b1_count = index.keys.select{|k| k.partition("~").last == "SF3B1"}
    iii sf3b1_count
  end
end

