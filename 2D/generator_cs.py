import argparse

parser = argparse.ArgumentParser(description='Array size')
parser.add_argument('--N', metavar='N', type=int,default = 4)
args = parser.parse_args()
COLS = vars(args)['N']
N = COLS
print(f"N={N}")
ROWS=COLS *2
LEN = (COLS + 1) * COLS / 2
LEN_S =  (COLS) * (COLS-1) / 2

def PE1_tail(id,oa,pb,cq,sq,fr):
    return(f"	PE1_tail(out_a[{oa}], pass_b[{pb}], out_cs_Q[{cq}], final_R[{fr}]);")

def PE1(id,fa,fb,oa,pc,ps,cq,sq,fr):
    if id ==0:
        return (f"	PE1({COLS}, Feed_A[0], Feed_B[0], out_a[0], PE1_pass_cs[0], out_cs_Q[0], final_R[0]);")
    else:
        return(f"	PE1({COLS-id}, out_a[{fa}], pass_b[{fb}], out_a[{oa}], PE1_pass_cs[{id}], out_cs_Q[{sq}], final_R[{fr}]);")

def PE2(idx,idy,fa,fb,pcin,psin,pa,pb,pcout,psout,fr):
    if idx == 0:
        if idy == 0:
            return (f"	PE2({COLS-idy-1}, Feed_A[{idy+1}], Feed_B[{idy+1}], PE1_pass_cs[{idx}], pass_a[{idy}], pass_b[{idy}], pass_cs[{idx}][{psout}], final_R[{fr}]);")
        else:
            return (f"	PE2({COLS-idy-1}, Feed_A[{idy+1}], Feed_B[{idy+1}], pass_cs[{idx}][{psin}], pass_a[{idy}], pass_b[{idy}], pass_cs[{idx}][{psout}], final_R[{fr}]);")
    else:
        if idy == 0:
            return(f"	PE2({COLS-idx-idy-1}, pass_a[{fa}], pass_b[{fb}], PE1_pass_cs[{idx}], pass_a[{pa}], pass_b[{pb}], pass_cs[{idx}][{psout}], final_R[{fr}]);")
        else:
            return(f"	PE2({COLS-idx-idy-1}, pass_a[{fa}], pass_b[{fb}], pass_cs[{idx}][{psin}], pass_a[{pa}], pass_b[{pb}], pass_cs[{idx}][{psout}], final_R[{fr}]);")

def PE2_tail(id,fa,fb,pcin,psin,pb,fr):
    if id ==0:
        return (f"	PE2_tail(Feed_A[{N-1}], Feed_B[{N-1}], pass_cs[{id}][{psin}], pass_b[{N-2}], final_R[{N-1}]);")
    if id == N-2:
        return (f"	PE2_tail(pass_a[{fa}], pass_b[{fb}], PE1_pass_cs[{id}], pass_b[{pb}], final_R[{fr}]);")
    return(f"	PE2_tail(pass_a[{fa}], pass_b[{fb}], pass_cs[{id}][{psin}], pass_b[{pb}], final_R[{fr}]);")

def PE3_head(cq,sq,ql,lqr,fql,fqr):
    return(f"	PE3_head(out_cs_Q[{sq}], Q_left[{ql}], Q_right[{lqr}], final_Q_left[{fql}], final_Q_right[{fqr}]);")

def PE3(id,cq,sq,ql,qr,oql,oqr,fql,fqr):
    return("	PE3({:},out_cs_Q[{:}],Q_left[{:}], Q_right[{:}],Q_left[{:}], Q_right[{:}], final_Q_left[{:}], final_Q_right[{:}]);".format(id,cq,ql,qr,oql,oqr,fql,fqr))
def PE3_tail(id,cq,sq,ql,qr,fql,fqr):
    return("	PE3_tail({:},out_cs_Q[{:}], Q_left[{:}], Q_right[{:}], final_Q_left[{:}], final_Q_right[{:}]);".format(id,cq,ql,qr,fql,fqr))


a = 1
pa = 0
pal = 0
par = N-2
pbl = 0
pbr = N-1
oa = 0
c = 0
d = 0;
fr = 0

PE1_list = []
PE2_list = []
PE3_list = []
pass_c = 0
pass_s = 0


for idx in range(0,N-1):
    pass_c = -1
    pass_s = -1
    PE2_list_y = []
    if idx ==1:
        pal = 0
        pbl = 0
        par = N-2
        pbr = N-1
    PE1_list.append(PE1(idx,idx-1,pbl,idx,pass_c,pass_s,idx,idx,fr))
    pbl+=1
    fr+=1
    for idy in range(N-idx-2):
        PE2_list_y.append(PE2(idx=idx,idy=idy,fa=pal,fb=pbl,pa=par,pb=pbr,fr=fr,psin = pass_c,pcin = pass_s,psout=pass_c+1,pcout=pass_s+1))
        pass_c += 1
        pass_s += 1
        pal+=1 
        par+=1
        pbl+=1
        pbr+=1
        fr+=1
    idy =0 if idx == COLS-2 else idy+1
    PE2_list_y.append(PE2_tail(idx,pal,pbl,pass_c,pass_s,pbr,fr))
    pass_c += 1
    pass_s += 1
    pal+=1
    pbl+=1
    pbr+=1
    fr+=1
    PE2_list.append(PE2_list_y)
PE1_list.append(PE1_tail(N-1,N-2,pbl,N-1,N-1,fr))

idx = 0
ql = 0
qr = 1
PE3_list.append(PE3_head(0,0,0,0,0,0))
for idx in range(1,N-1):
    PE3_list.append(PE3(idx,idx,idx,ql,ql,qr,qr,idx,idx))
    ql +=1
    qr +=1

PE3_list.append(PE3_tail(N-1,N-1,N-1,ql,ql,N-1,N-1))

def tick(idx):
    fw.write(PE1_list[idx]+'\n')
    for idy in range(idx+1):
        x = idx-idy-1
        y = 1+2*idy
        if x+y>N-2 or x<0:
            break
        fw.write(PE2_list[x][y]+'\n')

def tock(idx):
    fw.write(PE3_list[idx]+'\n')
    for idy in range(idx+1):
        x = idx-idy
        y = 0+2*idy
        if x+y>N-2 or x<0:
            break
        fw.write(PE2_list[x][y]+'\n')
        
def PE3():
    for pe3 in PE3_list:
        fw.write(pe3+'\n')
        
import os
copy = True
with open("dtqr2d.cpp","r") as f,open("temp","w") as fw:
    for line in f:
        if copy:
            fw.write(line)
        if "REPLACE START" in line:
            copy = False
            for i in range(N):
                tick(i)
                fw.write("\n")
                tock(i)
                fw.write("\n")
            # PE3()
        if "REPLACE END" in line:
            fw.write(line)
            copy = True
            
os.remove("dtqr2d.cpp")
os.rename("temp","dtqr2d.cpp")

with open("dtqr2d.h","r") as f,open("temp","w") as fw:
    for line in f:
        if "#define COLS" in line:
            line = f"#define COLS {N}\n"
        fw.write(line)
os.remove("dtqr2d.h")
os.rename("temp","dtqr2d.h")

                
print("finished")
test_matrix1 = "	MATRIX_T A1[LEN] = {"
test_matrix2 = "	MATRIX_T A2[LEN] = {"
LEN = int((1+N)*N/2)
import random

for idx in range(1,LEN):
    test_matrix1 = test_matrix1+str(random.randint(1,255))+" , "
    test_matrix2 = test_matrix2+str(random.randint(1,255))+" , "
test_matrix1 = test_matrix1 + str(random.randint(1,255))+"};"
test_matrix2 = test_matrix2 + str(random.randint(1,255))+"};"
print(test_matrix1)
print(test_matrix2)
