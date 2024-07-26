# Usage: python3 densities.py <my OBD filename> <my TBD filename> <trdens filename>

import sys
import math

# Convert number to Fortran format
def convert(number):
    a = str('{:.7E}'.format(number))
    if a[len(a)-3:] == '-01':
        sign = '+'
    else:
        sign = a[len(a)-3]
    if number == 0.0:
        end = str(int(a[len(a)-2:]))
    elif a[len(a)-3] == '+':
        end = str(int(a[len(a)-2:])+1)
    else:
        end = str(int(a[len(a)-2:])-1)
    if len(end) == 1:
        end = '0' + end
    if a[0] == '-':
        converted = '-0.' + a[1] + a[3:11] + sign + end
    else:
        converted = ' 0.' + a[0] + a[2:10] + sign + end
    return converted

# Read my OBDRMEs
obd_file = open(sys.argv[1], 'r')
obdrmes={} # key is (Np,lp,jjp,N,l,jj,J0), value is (proton OBDRME, neutron OBDRME)
line_count = 0
for line in obd_file:
    line_count += 1
    i = 0
    f = 0
    k = 0
    numbers = []
    for char in line:
        f += 1
        if char == ' ' or f == len(line):
            k += 1
            if k <= 7:
                numbers.append(int(line[i:f]))
            else:
                numbers.append(float(line[i:f]))
            i=f
    Np = numbers[0]
    lp = numbers[1]
    jjp = numbers[2]
    N = numbers[3]
    l = numbers[4]
    jj = numbers[5]
    J0 = numbers[6]
    p_obd = numbers[7]
    n_obd = numbers[8]
    key = (Np,lp,jjp,N,l,jj,J0)
    obdrmes[key] = (p_obd,n_obd)
obd_file.close()

# Read my TBDRMEs
tbd_file = open(sys.argv[2], 'r')
tbdrmes = {} # key is (N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0), value is (proton TBDRME, neutron TBDRME, pn TBDRME)
line_count = 0
for line in tbd_file:
    line_count += 1
    i = 0
    f = 0
    k = 0
    numbers = []
    for char in line:
        f += 1
        if char == ' ' or f == len(line):
            k += 1
            if k <= 15:
                numbers.append(int(line[i:f]))
            else:
                numbers.append(float(line[i:f]))
            i=f
    N1 = numbers[0]
    l1 = numbers[1]
    jj1 = numbers[2]
    N2 = numbers[3]
    l2 = numbers[4]
    jj2 = numbers[5]
    N3 = numbers[6]
    l3 = numbers[7]
    jj3 = numbers[8]
    N4 = numbers[9]
    l4 = numbers[10]
    jj4 = numbers[11]
    Jf = int(numbers[12]/2)
    Ji = int(numbers[13]/2)
    J0 = int(numbers[14]/2)
    p_tbd = numbers[15]
    n_tbd = numbers[16]
    pn_tbd = numbers[17]
    key = (N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0)
    if N1 == N2 and l1 == l2 and jj1 == jj2:
        factor = math.sqrt(2.0)
    else:
        factor = 1.0
    if N3 == N4 and l3 == l4 and jj3 == jj4:
        factor*=math.sqrt(2.0)
    tbdrmes[key]=(numbers[15]/factor,numbers[16]/factor,numbers[17]/factor)
tbd_file.close()

# Replace OBDRMEs and TBDRMEs in trdens file with my RMEs
trdens_file = open(sys.argv[3], 'r')
sp_states = []
tb_labels_by_index = {} # key is index, value is (N1,l1,jj1,N2,l2,jj2,J,T)
for line in trdens_file:
    if len(line) == 29:
        if line[8]=='n':
            n = int(line[10:13])
            l = int(line[17:20])
            jj = int(line[24:26])
            N = 2*n + l
            sp_states.append((N,l,jj))
        else:
            index = int(line[2:7])
            a = int(line[17:20])
            b = int(line[20:23])
            J = int(line[23:26])
            T = int(line[27])
            N1 = sp_states[a-1][0]
            l1 = sp_states[a-1][1]
            jj1 = sp_states[a-1][2]
            N2 = sp_states[b-1][0]
            l2 = sp_states[b-1][1]
            jj2 = sp_states[b-1][2]
            tb_labels_by_index[index]=(N1,l1,jj1,N2,l2,jj2,J,T)
        print(line[0:len(line)-1])
    elif len(line) == 12:
        J0 = int(line[8:11])
        print(line[0:len(line)-1])
    elif len(line) == 42:
        a = int(line[0:4])
        b = int(line[5:9])
        Np = sp_states[a-1][0]
        lp = sp_states[a-1][1]
        jjp = sp_states[a-1][2]
        N = sp_states[b-1][0]
        l = sp_states[b-1][1]
        jj = sp_states[b-1][2]
        p_obd = obdrmes[(Np,lp,jjp,N,l,jj,J0)][0]
        n_obd = obdrmes[(Np,lp,jjp,N,l,jj,J0)][1]
        np = (Np-lp)/2
        n = (N-l)/2
        if((N+np+n)%2 != 0):
            p_obd = -p_obd
            n_obd = -n_obd
        proton_obdrme = float(line[10:26])
        neutron_obdrme = float(line[26:])
#        if(p_obd*proton_obdrme<0.0):
#            p_obd = -p_obd
#        if(n_obd*neutron_obdrme<0.0):
#            n_obd = -n_obd
#        if(abs(p_obd-proton_obdrme)>1.0e-5):
#            print('ERROR: proton OBDRMEs:',p_obd,proton_obdrme)
#        if(abs(n_obd-neutron_obdrme)>1.0e-5):
#            print('ERROR: neutron OBDRMEs:',n_obd,neutron_obdrme)
        print(line[0:10]+convert(p_obd)+' '+convert(n_obd))
    elif line[0:2] == 'tb':
        index12 = int(line[3:8])
        index34 = int(line[9:15])
        pn_tbdrme = float(line[16:31])
        proton_tbdrme = float(line[32:47])
        neutron_tbdrme = float(line[48:])
        N1 = tb_labels_by_index[index12][0]
        l1 = tb_labels_by_index[index12][1]
        jj1 = tb_labels_by_index[index12][2]
        N2 = tb_labels_by_index[index12][3]
        l2 = tb_labels_by_index[index12][4]
        jj2 = tb_labels_by_index[index12][5]
        Jf = tb_labels_by_index[index12][6]
        Tf = tb_labels_by_index[index12][7]
        N3 = tb_labels_by_index[index34][0]
        l3 = tb_labels_by_index[index34][1]
        jj3 = tb_labels_by_index[index34][2]
        N4 = tb_labels_by_index[index34][3]
        l4 = tb_labels_by_index[index34][4]
        jj4 = tb_labels_by_index[index34][5]
        Ji = tb_labels_by_index[index34][6]
        Ti = tb_labels_by_index[index34][7]
        if Tf == 1 and Ti == 1:
            p_tbd = tbdrmes[(N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0)][0]
            n_tbd = tbdrmes[(N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0)][1]
            pn_tbd = tbdrmes[(N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0)][2]
            if((Ji-(jj3+jj4)/2)%2 == 0):
                pn_tbd = pn_tbd-tbdrmes[(N1,l1,jj1,N2,l2,jj2,N4,l4,jj4,N3,l3,jj3,Jf,Ji,J0)][2]
            else:
                pn_tbd = pn_tbd+tbdrmes[(N1,l1,jj1,N2,l2,jj2,N4,l4,jj4,N3,l3,jj3,Jf,Ji,J0)][2]
            if((Jf-(jj1+jj2)/2)%2 == 0):
                pn_tbd = pn_tbd-tbdrmes[(N2,l2,jj2,N1,l1,jj1,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0)][2]
            else:
                pn_tbd = pn_tbd+tbdrmes[(N2,l2,jj2,N1,l1,jj1,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0)][2]
            if((Ji+Jf-(jj3+jj4+jj1+jj2)/2)%2 == 0):
                pn_tbd = pn_tbd+tbdrmes[(N2,l2,jj2,N1,l1,jj1,N4,l4,jj4,N3,l3,jj3,Jf,Ji,J0)][2]
            else:
                pn_tbd = pn_tbd-tbdrmes[(N2,l2,jj2,N1,l1,jj1,N4,l4,jj4,N3,l3,jj3,Jf,Ji,J0)][2]
            pn_tbd = pn_tbd/2
        elif Tf == 0 and Ti == 0:
            p_tbd = 0.0
            n_tbd = 0.0
            pn_tbd = -tbdrmes[(N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0)][2]
            if((Ji-(jj3+jj4)/2)%2 == 0):
                pn_tbd = pn_tbd-tbdrmes[(N1,l1,jj1,N2,l2,jj2,N4,l4,jj4,N3,l3,jj3,Jf,Ji,J0)][2]
            else:
                pn_tbd = pn_tbd+tbdrmes[(N1,l1,jj1,N2,l2,jj2,N4,l4,jj4,N3,l3,jj3,Jf,Ji,J0)][2]
            if((Jf-(jj1+jj2)/2)%2 == 0):
                pn_tbd = pn_tbd-tbdrmes[(N2,l2,jj2,N1,l1,jj1,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0)][2]
            else:
                pn_tbd = pn_tbd+tbdrmes[(N2,l2,jj2,N1,l1,jj1,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0)][2]
            if((Ji+Jf-(jj3+jj4+jj1+jj2)/2)%2 == 0):
                pn_tbd = pn_tbd-tbdrmes[(N2,l2,jj2,N1,l1,jj1,N4,l4,jj4,N3,l3,jj3,Jf,Ji,J0)][2]
            else:
                pn_tbd = pn_tbd+tbdrmes[(N2,l2,jj2,N1,l1,jj1,N4,l4,jj4,N3,l3,jj3,Jf,Ji,J0)][2]
            pn_tbd = pn_tbd/2
        elif Tf == 0 and Ti == 1:
            p_tbd = 0.0
            n_tbd = 0.0
            pn_tbd = tbdrmes[(N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0)][2]
            if((Ji-(jj3+jj4)/2)%2 == 0):
                pn_tbd = pn_tbd-tbdrmes[(N1,l1,jj1,N2,l2,jj2,N4,l4,jj4,N3,l3,jj3,Jf,Ji,J0)][2]
            else:
                pn_tbd = pn_tbd+tbdrmes[(N1,l1,jj1,N2,l2,jj2,N4,l4,jj4,N3,l3,jj3,Jf,Ji,J0)][2]
            if((Jf-(jj1+jj2)/2)%2 == 0):
                pn_tbd = pn_tbd+tbdrmes[(N2,l2,jj2,N1,l1,jj1,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0)][2]
            else:
                pn_tbd = pn_tbd-tbdrmes[(N2,l2,jj2,N1,l1,jj1,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0)][2]
            if((Ji+Jf-(jj3+jj4+jj1+jj2)/2)%2 == 0):
                pn_tbd = pn_tbd-tbdrmes[(N2,l2,jj2,N1,l1,jj1,N4,l4,jj4,N3,l3,jj3,Jf,Ji,J0)][2]
            else:
                pn_tbd = pn_tbd+tbdrmes[(N2,l2,jj2,N1,l1,jj1,N4,l4,jj4,N3,l3,jj3,Jf,Ji,J0)][2]
            pn_tbd = pn_tbd/2
        else:
            p_tbd = 0.0
            n_tbd = 0.0
            pn_tbd = -tbdrmes[(N1,l1,jj1,N2,l2,jj2,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0)][2]
            if((Ji-(jj3+jj4)/2)%2 == 0):
                pn_tbd = pn_tbd-tbdrmes[(N1,l1,jj1,N2,l2,jj2,N4,l4,jj4,N3,l3,jj3,Jf,Ji,J0)][2]
            else:
                pn_tbd = pn_tbd+tbdrmes[(N1,l1,jj1,N2,l2,jj2,N4,l4,jj4,N3,l3,jj3,Jf,Ji,J0)][2]
            if((Jf-(jj1+jj2)/2)%2 == 0):
                pn_tbd = pn_tbd+tbdrmes[(N2,l2,jj2,N1,l1,jj1,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0)][2]
            else:
                pn_tbd = pn_tbd-tbdrmes[(N2,l2,jj2,N1,l1,jj1,N3,l3,jj3,N4,l4,jj4,Jf,Ji,J0)][2]
            if((Ji+Jf-(jj3+jj4+jj1+jj2)/2)%2 == 0):
                pn_tbd = pn_tbd+tbdrmes[(N2,l2,jj2,N1,l1,jj1,N4,l4,jj4,N3,l3,jj3,Jf,Ji,J0)][2]
            else:
                pn_tbd = pn_tbd-tbdrmes[(N2,l2,jj2,N1,l1,jj1,N4,l4,jj4,N3,l3,jj3,Jf,Ji,J0)][2]
            pn_tbd = pn_tbd/2
        n1 = (N1-l1)/2
        n2 = (N2-l2)/2
        n3 = (N3-l3)/2
        n4 = (N4-l4)/2
        if((N3+N4+n1+n2+n3+n4+1)%2 != 0):
            p_tbd = -p_tbd
            n_tbd = -n_tbd
            pn_tbd = -pn_tbd
        pn_tbdrme = float(line[16:31])
        proton_tbdrme = float(line[32:47])
        neutron_tbdrme = float(line[48:])
#        if(p_tbd*proton_tbdrme<0.0):
#            p_tbd = -p_tbd
#        if(n_tbd*neutron_tbdrme<0.0):
#            n_tbd = -n_tbd
#        if(pn_tbd*pn_tbdrme<0.0):
#            pn_tbd = -pn_tbd
#        if(abs(p_tbd-proton_tbdrme)>1.0e-5):
#            print('ERROR: proton TBDRMEs:',p_tbd,proton_tbdrme)
#        if(abs(n_tbd-neutron_tbdrme)>1.0e-5):
#            print('ERROR: neutron TBDRMEs:',n_tbd,neutron_tbdrme)
#        if(abs(pn_tbd-pn_tbdrme)>1.0e-5):
#            print('ERROR: proton-neutron TBDRMEs:',pn_tbd,pn_tbdrme)
        print(line[0:16]+convert(pn_tbd)+' '+convert(p_tbd)+' '+convert(n_tbd))
    else:
        print(line[0:len(line)-1])
trdens_file.close()
