import os
import sys

blas_template = \
"""
#ifndef OP2_DG_CUDA
#if DG_DOUBLE == 1
cblas_dgemv({},{},{},{},{},{},{},{}, 1,{},{}, 1);
#else
cblas_sgemv({},{},{},{},{},{},{},{}, 1,{},{}, 1);
#endif
#else
for(int i = 0; i < {}; i++) {{
  ({})[i] *= {};
  for(int j = 0; j < {}; j++) {{
    int ind = DG_MAT_IND({},{},{},{});
    ({})[i] += {} * ({})[ind] * ({})[j];
  }}
}}
#endif
"""

def replace_blas_kernels(input_str):
    out_str = input_str
    index = out_str.find("op2_in_kernel_gemv")
    col_maj = True
    while index != -1:
        end_ind = out_str.find(")", index)
        args_str = out_str[out_str.find("(", index) + 1 : end_ind]
        args_str = args_str.split(",")
        col_row_maj_str = "CblasColMajor"
        if not col_maj:
            col_row_maj_str = "CblasRowMajor"
        transpose_str = " CblasNoTrans"
        ind0 = "i"
        ind1 = "j"
        if "true" in args_str[0]:
            transpose_str = " CblasTrans"
            ind0 = "j"
            ind1 = "i"
        blas_call = blas_template.format(col_row_maj_str, transpose_str, args_str[1], \
            args_str[2], args_str[3], args_str[4], args_str[5], args_str[6], \
            args_str[7], args_str[8], col_row_maj_str, transpose_str, args_str[1], \
            args_str[2], args_str[3], args_str[4], args_str[5], args_str[6], \
            args_str[7], args_str[8], args_str[1], args_str[8], args_str[7], \
            args_str[2], ind0, ind1, args_str[1], args_str[2], args_str[8], args_str[3], \
            args_str[4], args_str[6])
        out_str = out_str[0 : index] + blas_call + out_str[end_ind + 2 :]
        index = out_str.find("op2_in_kernel_gemv")
    return out_str

dim = sys.argv[1]
order = sys.argv[2]

dg_np        = ""
dg_npf       = ""
dg_cub_np    = ""
dg_g_np      = ""
dg_gf_np     = ""
ls_sample_np = ""
dg_num_faces = ""
dg_order     = order
dg_num_constants = ""
fp_type = "d"

# Get DG order from command line args
if dim == "2":
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
        sys.exit(-1)
elif dim == "3":
    order_int = int(order)
    dg_np     = str(int((order_int + 1) * (order_int + 2) * (order_int + 3) / 6))
    dg_npf    = str(int((order_int + 1) * (order_int + 2) / 2))
    dg_num_faces = "4"
    dg_num_constants = "2"
    ls_sample_np = "13"

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
    newdata = replace_blas_kernels(newdata)
    if "CMakeLists" not in f:
        if fp_type == "d":
            newdata = newdata.replace("DG_FP_STR", "\"double\"")
            newdata = newdata.replace("DG_FP", "double")
            newdata = newdata.replace("DG_MPI_FP", "MPI_DOUBLE")
        else:
            newdata = newdata.replace("DG_FP_STR", "\"float\"")
            newdata = newdata.replace("DG_FP", "float")
            newdata = newdata.replace("DG_MPI_FP", "MPI_FLOAT")
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
