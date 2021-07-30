import os

# Order 4
dg_np     = "15"
dg_npf    = "5"
dg_cub_np = "46"
dg_g_np   = "21"
dg_gf_np  = "7"
# Order 3
# dg_np     = "10"
# dg_npf    = "4"
# dg_cub_np = "36"
# dg_g_np   = "18"
# dg_gf_np  = "6"
# Order 2
# dg_np     = "6"
# dg_npf    = "3"
# dg_cub_np = "16"
# dg_g_np   = "12"
# dg_gf_np  = "4"
# Order 1
# dg_np     = "3"
# dg_npf    = "2"
# dg_cub_np = "12"
# dg_g_np   = "9"
# dg_gf_np  = "3"

inputfiles = []

for dirpath, _, filenames in os.walk("src"):
    for f in filenames:
        tmp  = dirpath + "/" + f
        tmp2 = tmp.split("/")
        tmp3 = "/".join(tmp2[1:])
        inputfiles.append(tmp3)

for f in inputfiles:
    filedata = None
    with open("src/" + f, "r") as file:
        filedata = file.read()

    newdata = filedata.replace("DG_NPF", dg_npf)
    newdata = newdata.replace("DG_NP", dg_np)
    newdata = newdata.replace("DG_CUB_NP", dg_cub_np)
    newdata = newdata.replace("DG_G_NP", dg_g_np)
    newdata = newdata.replace("DG_GF_NP", dg_gf_np)

    with open("gen/" + f, "w") as file:
        file.write(newdata)
