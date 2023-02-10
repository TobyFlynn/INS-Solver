import os
import sys

dim = sys.argv[1]
order = sys.argv[2]

dg_np        = ""
dg_npf       = ""
dg_cub_np    = ""
dg_g_np      = ""
dg_gf_np     = ""
ls_sample_np = ""
dg_num_faces = ""
dg_num_constants = ""
dg_order     = sys.argv[1]

# Get DG order from command line args
if dim == "2":
    dg_order  = order
    order_int = int(order)
    dg_np     = str(int((order_int + 1) * (order_int + 2) / 2))
    dg_npf    = str(order_int + 1)
    ls_sample_np = "10"
    dg_num_faces = "3"
    dg_num_constants = "5"
    if order == "1":
        dg_cub_np    = "12"
        dg_g_np      = "9"
        dg_gf_np     = "3"
    elif order == "2":
        dg_cub_np    = "16"
        dg_g_np      = "12"
        dg_gf_np     = "4"
    elif order == "3":
        dg_cub_np    = "36"
        dg_g_np      = "18"
        dg_gf_np     = "6"
    elif order == "4":
        dg_cub_np    = "46"
        dg_g_np      = "21"
        dg_gf_np     = "7"
    else:
        print("This order of DG is not supported yet...exiting")
        sys.exit()
elif dim == "3":
    dg_order  = order
    order_int = int(order)
    dg_np     = str(int((order_int + 1) * (order_int + 2) * (order_int + 3) / 6))
    dg_npf    = str(int((order_int + 1) * (order_int + 2) / 2))
    dg_num_faces = "4"
    dg_num_constants = "2"

inputfiles = []

for dirpath, _, filenames in os.walk("src"):
    for f in filenames:
        if f[0] != '.':
            tmp  = dirpath + "/" + f
            tmp2 = tmp.split("/")
            tmp3 = "/".join(tmp2[1:])
            inputfiles.append(tmp3)

for f in inputfiles:
    filedata = None
    with open("src/" + f, "r") as file:
        filedata = file.read()

    newdata = filedata
    if "CMakeLists" not in f:
        newdata = newdata.replace("DG_NPF", dg_npf)
        newdata = newdata.replace("DG_NP", dg_np)
        newdata = newdata.replace("DG_CUB_NP", dg_cub_np)
        newdata = newdata.replace("DG_G_NP", dg_g_np)
        newdata = newdata.replace("DG_GF_NP", dg_gf_np)
        newdata = newdata.replace("LS_SAMPLE_NP", ls_sample_np)
        newdata = newdata.replace("DG_NUM_FACES", dg_num_faces)
        newdata = newdata.replace("DG_NUM_CONSTANTS", dg_num_constants)
        newdata = newdata.replace("DG_ORDER", dg_order)

    if dim == "2":
        with open("gen2d/" + f, "w") as file:
            file.write(newdata)
    elif dim == "3":
        with open("gen3d/" + f, "w") as file:
            file.write(newdata)
