# edit the base class file name
cf_basename = "../../CMS3" # .cc, .h
# copy the pattern for the following branch (make sure the type is correct)
branch_to_copy = "els_miniIso_uncor"
branch_to_copy_type = "vector<float>"
# and use this pattern to squeeze in the following branches (and associated types)
branches_to_add = [
        # ("mus_miniRelIso_chg", "vector<float>"),
        # ("mus_miniRelIso_all", "vector<float>"),
        # ("els_miniRelIso_chg", "vector<float>"),
        # ("els_miniRelIso_all", "vector<float>"),
        # ("els_ptRatio", "vector<float>"),
        # ("els_ptRel", "vector<float>"),
        # ("els_jetBTagCSV", "vector<float>"),
        # ("mus_ptRatio", "vector<float>"),
        # ("mus_ptRel", "vector<float>"),
        # ("mus_jetBTagCSV", "vector<float>"),
        # ("els_jetNDauChargedMVASel", "vector<int>"),
        # ("mus_jetNDauChargedMVASel", "vector<int>"),

        ("genHEPMCweight", "float"),
        ("genHEPMCweight_2016", "float"),
        ("gen_LHEweight_QCDscale_muR1_muF1", "float"),
        ("gen_LHEweight_QCDscale_muR1_muF2", "float"),
        ("gen_LHEweight_QCDscale_muR1_muF0p5", "float"),
        ("gen_LHEweight_QCDscale_muR2_muF1", "float"),
        ("gen_LHEweight_QCDscale_muR2_muF2", "float"),
        ("gen_LHEweight_QCDscale_muR2_muF0p5", "float"),
        ("gen_LHEweight_QCDscale_muR0p5_muF1", "float"),
        ("gen_LHEweight_QCDscale_muR0p5_muF2", "float"),
        ("gen_LHEweight_QCDscale_muR0p5_muF0p5", "float"),
        ("gen_LHEweight_PDFVariation_Up", "float"),
        ("gen_LHEweight_PDFVariation_Dn", "float"),
        ("gen_LHEweight_AsMZ_Up", "float"),
        ("gen_LHEweight_AsMZ_Dn", "float"),
        ("gen_LHEweight_PDFVariation_Up_2016", "float"),
        ("gen_LHEweight_PDFVariation_Dn_2016", "float"),
        ("gen_LHEweight_AsMZ_Up_2016", "float"),
        ("gen_LHEweight_AsMZ_Dn_2016", "float"),

        ]

with open("{}_new.h".format(cf_basename), "w") as fhout:
    with open("{}.h".format(cf_basename), "r") as fhin:
        for line in fhin:
            if "bool {}_isLoaded;".format(branch_to_copy) in line:
                fhout.write(line)
                for bname,btype in branches_to_add:
                    fhout.write("  {} {}_;\n".format(btype,bname))
                    fhout.write("  TBranch *{}_branch;\n".format(bname))
                    fhout.write("  bool {}_isLoaded;\n".format(bname))
            elif "const {} &{}();".format(branch_to_copy_type,branch_to_copy) in line:
                fhout.write(line)
                for bname,btype in branches_to_add:
                    fhout.write("  const {} &{}();\n".format(btype,bname))
            else:
                fhout.write(line)


with open("{}_new.cc".format(cf_basename), "w") as fhout:
    with open("{}.cc".format(cf_basename), "r") as fhin:
        for line in fhin:
            if "{}_branch = 0;".format(branch_to_copy) == line.strip():
                for bname,btype in branches_to_add:
                    fhout.write("""
  {0}_branch = 0;
  if (tree->GetAlias(\"{0}\") != 0) {{
    {0}_branch = tree->GetBranch(tree->GetAlias(\"{0}\"));
    if ({0}_branch) {{ {0}_branch->SetAddress(&{0}_); }}
  }}
""".format(bname))
                fhout.write(line)

            elif "{}_isLoaded = false;".format(branch_to_copy) == line.strip():
                fhout.write(line)
                for bname,btype in branches_to_add:
                    fhout.write("  {}_isLoaded = false;\n".format(bname))

            elif "const {0} &{1}() {{ return cms3.{1}(); }}".format(branch_to_copy_type,branch_to_copy) == line.strip():
                fhout.write(line)
                for bname,btype in branches_to_add:
                    fhout.write("  const {0} &{1}() {{ return cms3.{1}(); }}\n".format(btype,bname))

            elif "if ({0}_branch != 0) {0}();".format(branch_to_copy) == line.strip():
                fhout.write(line)
                for bname,btype in branches_to_add:
                    fhout.write("  if ({0}_branch != 0) {0}();\n".format(bname))

            elif "const {} &CMS3::{}() {{".format(branch_to_copy_type,branch_to_copy) == line.strip():
                for bname,btype in branches_to_add:
                    fhout.write("""
const {0} &CMS3::{1}() {{
  if (not {1}_isLoaded) {{
    if ({1}_branch != 0) {{
      {1}_branch->GetEntry(index);
    }} else {{
      printf("branch {1}_branch does not exist!\\n");
      exit(1);
    }}
    {1}_isLoaded = true;
  }}
  return {1}_;
}}
""".format(btype,bname))
                fhout.write(line)

            else:
                fhout.write(line)

