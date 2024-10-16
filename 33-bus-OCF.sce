//----------------------------------------------
// IEEE 34-Bus
// @autor: Adolfo Blengini Neto
// @autor: Marcius Fabius Henriques source Carvalho
// @autor: Secundino Soares Filho
// OCF model
// @power distribution system: IEEE 34-bu, DERs = bus 18-28, Voltage Control, Line Limits
//----------------------------------------------
clear; clc;

data_bus = [
//  Nr type VMr   VMi    PC      Qc        Pg    Qg  BshT  Qmax  Qmin    Pmax   Pmin      Zr             Zi           Ir          Ii        Yr       Yi
//          [pu]  [pu] [MW]     [Mvar]    [MW] [Mvar]    [Mvar] [Mvar]  [MW]   [MW]
    1   3   1.00    0    0         0         0    0    0    4.0   -4.0     4.1   0        9.60923      18.58736    0.08690     -0.04493      0        0 
    2   0   1.00    0    0         0         0    0    0    0       0      0     0        0            0           0            0            0        0
    3   0   1.00    0    0         0         0    0    0    0       0      0     0       18.89169      36.54081    0            0            0        0
    4   0   1.00    0    0         0         0    0    0    0       0      0     0        0            0           0            0            0        0
    5   0   1.00    0    0         0         0    0    0    0       0      0     0        0            0           0            0            0        0
    6   0   1.00    0    0         0         0    0    0    0       0      0     0        0            0           0            0            0        0
    7   0   1.00    0    0         0         0    0    0    0       0      0     0        0            0           0            0            0        0
    8   0   1.00    0    0         0         0    0    0    0       0      0     0        769.2308     1500.0      0            0            0        0
    9   0   1.00    0    0.1130   -0.23357   0    0    0    0       0      0     0        0            0           0            0            0        0
    10  0   1.00    0    0.45177  -0.23357   0    0    0    0       0      0     0        0            0           0            0            0        0
    11  0   1.00    0    0         0         0    0    0    0       0      0     0        0            0           0            0            0        0
    12  0   1.00    0    0         0         0    0    0    0       0      0     0        7.15478      13.83764    0.00927     -0.00480      0        0
    13  0   1.00    0    0         0         0    0    0    0       0      0     0        0            0           0            0            0        0
    14  0   1.00    0    0.02060  -0.01067   0    0    0    0       0      0     0        0            0           0            0            0        0
    15  0   1.00    0    0         0         0    0    0    0       0      0     0        0            0           0            0            0        0
    16  0   1.00    0    0         0         0    0    0    0       0      0     0        80.86253     156.25000   0            0            0        0
    17  0   1.00    0    0         0         0    0    0    0       0      0     0        0            0           0            0            0        0
//  18  0   1.00    0    0         0         0    0    0    0       0      0     0        0            0           0            0            0        0
    18  2   1.00    0    0         0         0    0    0    0       0      0.5   0.1      0            0           0            0            0        0
    19  0   1.00    0    0.02227  -0.01150   0    0    0    0       0      0     0        277.77780     535.71430  0.01783     -0.00923      0        0
    20  0   1.00    0    0         0         0    0    0    0       0      0     0        0            0           0            0            0        0
    21  0   1.00    0    0.27000  -0.21620   0    0    0    0       0      0     0        0            0           0            0            0        0
    22  0   1.00    0    0.01540  -0.00797   0    0    0    0       0      0     0       23.90438      46.22496    0.04273     -0.02210      0        0
    23  0   1.00    0    0         0         0    0    0    0       0      0     0        0            0           0            0            0        0
    24  0   1.00    0    0.05220  -0.02697   0    0    0    0       0      0     0       14.38159      27.82932    0.37050     -0.19153      0        0
    25  0   1.00    0    0.03040  -0.01570   0    0    0    0       0      0     0        0            0           0            0            0        0
    26  0   1.00    0    0.33440  -0.06830   0    0    0    0       0      0     0       12.20008      23.60346    0.07410     -0.03830      0        1.0
    27  0   1.00    0    0         0         0    0    0    0       0      0     0       13.26260      25.64103    0            0            0        0
//  28  0   1.00    0    0.19450  -0.15570   0    0    0    0       0      0     0        0            0           0            0            0        1.5
    28  2   0.98    0    0.19450   -0.15570  0    0    0    0.8     0     2.0    0        0            0           0            0            0        1.5
    29  0   1.00    0    0.29013  -0.20657   0    0    0    0       0      0     0       28.43602      55.04587    0.14017     -0.07247      0        0
    30  0   1.00    0    0.05830  -0.03013   0    0    0    0       0      0     0       13.75516      26.61934    0            0            0        0
    31  0   1.00    0    0.08860  -0.07090   0    0    0    0       0      0     0        0            0           0            0            0        0
    32  0   1.00    0    0.09203  -0.04757   0    0    0    0       0      0     0        0            0           0            0            0        0
    33  0   1.00    0    0.        0.        0    0    0    0       0      0     0        0            0           0            0            0        0
];

//type 3: Vth. slack
//type 2: PV
//type 0: PQ

data_line = [
// source sourcest.  R          X          Bsh    Tap     Fi
//              [pu]       [pu]       [pu]
    1     2     0.00077    0.00051    0.0    0.0    0.0 
    2     3     0.00966    0.00644    0.0    0.0    0.0 
    3     4     0.00175    0.00116    0.0    0.0    0.0 
    3     5     0.01125    0.00750    0.0    0.0    0.0 
    5     6     0.00891    0.00594    0.0    0.0    0.0 
    6     7     0.00003    0.00002    0.0    0.0    0.0 
    7     8     0.00009    0.00006    0.0    0.0    0.0 
    8     9     0.00051    0.00034    0.0    0.0    0.0 
    9     10    0.01444    0.00963    0.0    0.0    0.0 
    10    11    0.01444    0.00963    0.0    0.0    0.0 
    8     12    0.00306    0.00204    0.0    0.0    0.0 
    12    13    0.00306    0.00204    0.0    0.0    0.0 
    12    14    0.00025    0.00016    0.0    0.0    0.0 
    14    15    0.00613    0.00408    0.0    0.0    0.0 
    15    16    0.00015    0.00010    0.0    0.0    0.0 
    16    17    0.00015    0.00010    0.0    0.0    0.0 
    16    18    0.01104    0.00736    0.0    0.0    0.0 
    18    19    0.00003    0.00002    0.0    0.0    0.0 
    19    20    0.00003    0.00002    0.0    0.0    0.0 
    20    21    0.00316    0.00211    0.0    0.0    0.0 
    19    22    0.00147    0.00098    0.0    0.0    0.0 
    22    23    0.00048    0.00032    0.0    0.0    0.0 
    22    24    0.00174    0.00116    0.0    0.0    0.0 
    24    25    0.00008    0.00005    0.0    0.0    0.0 
    25    26    0.00040    0.00027    0.0    0.0    0.0 
    26    27    0.00109    0.00072    0.0    0.0    0.0 
    27    28    0.00015    0.00010    0.0    0.0    0.0 
    24    29    0.00060    0.00040    0.0    0.0    0.0 
    29    30    0.00080    0.00053    0.0    0.0    0.0 
    30    31    0.00025    0.00017    0.0    0.0    0.0 
    30    32    0.00008    0.00005    0.0    0.0    0.0 
    32    33    0.00145    0.00097    0.0    0.0    0.0 
];
 
Sbase = 1;//Base 1
Vbase = 132;

vub = 100
vlb = -100


// -------------------------------
//  STEP1 : Read the system input
// ------------------------------

// read data bus
num = data_bus(:,1);
tipo = data_bus(:,2);
VM = data_bus(:,3);
Th = data_bus(:,4);
Pc = data_bus(:,5);
Qc = data_bus(:,6);
Pg = data_bus(:,7);
Qg = data_bus(:,8);
bsht = data_bus(:,9);
Qmax = data_bus(:,10);
Qmin = data_bus(:,11);
Pmax = data_bus(:,12);
Pmin = data_bus(:,13);
Zr = data_bus(:,14);
Zi = data_bus(:,15);
Ir = data_bus(:,16);
Ii = data_bus(:,17);
Yr = data_bus(:,18);
Yi = data_bus(:,19);


//read data line
source = data_line(:,1);
destination = data_line(:,2);
r = data_line(:,3);
x = data_line(:,4);

// system dimensions
nb = length(num);
nr = length(source);

// P and Q
Pesp = zeros(nb,1);
Qesp = zeros(nb,1);
Pcalc = zeros(nb,1);
Qcalc = zeros(nb,1);

// convert to P.U
for i = 1:nb
    Pc(i) = Pc(i)/Sbase;
    Qc(i) = Qc(i)/Sbase;
    Pg(i) = Pg(i)/Sbase;
    Qg(i) = Qg(i)/Sbase;
    Zr(i) = Zr(i)/Sbase;
    Zi(i) = Zi(i)/Sbase;
    Ir(i) = Ir(i)/Sbase;
    Ii(i) = Ii(i)/Sbase;
    Yr(i) = Yr(i)/Sbase;
    Yi(i) = Yi(i)/Sbase;
    
    Qmax(i) = Qmax(i)/Sbase;
    Qmin(i) = Qmin(i)/Sbase;
    Pmax(i) = Pmax(i)/Sbase;
    Pmin(i) = Pmin(i)/Sbase;
    Pesp(i) = Pg(i) - Pc(i);
    Qesp(i) = Qg(i) - Qc(i);
end

// ------------------------------
//  STEP2 : Create Variables
// ------------------------------

// matriz A
A = zeros(2*nb,4*nr);

// matriz V
V = zeros(2*nb-2,4*nr);

// matriz w
w = zeros(4 * nr,1);

//bus complex voltage vector
Vr_inicial = ones(nb,1);
Vi_inicial = zeros(nb,1);
Vr_inicial(1) = 0; 

//bus net complex current vector(nb x 1).
Ibus = zeros(nb,1);

// complex voltage at reference bus
Vf = 1.0 + %i * 0.0;

// number of sourceRs
gd = 0;

// the convergence factor 
erro = 1e-4;



e = zeros(nr,1);
meshdestination = zeros(nr,1);
for l = 1:nr
    m = destination(l);
    for (j = 1: l)
        if (m == e(j))
            disp ("****** Network is mesh *********");
            meshdestination(j) = m;
        end
    end
    e(l) = m;
end


// -----------------------------
//  STEP3 : Create A Matrix
// -----------------------------
//  A Matrix line = source-sourcestination  
//  -
//  -1 (source)
//   1 (sourcestination)
// - 
for l = 1:nr
    k = source(l);
    m = destination(l);
    
    //Real
    A(k,l) = -1
    A(m,l) =  1
    
    //Imag.
    A(k+nb,l+nr) = -1
    A(m+nb,l+nr) =  1
    
    //costs
    w(l) = r(l);
    w(l+nr) = r(l);
end


// ----------------------------
// STEP 4 : Create External / Internal Loops 
// ---------------------------
// - - - - - - - - -
// R    -X    1
// X     R         1
// - - - - - - - - -
lei = 1;
for t = 1:nr
    
    //line by line
    k = source(t);
    m = destination(t);
    
    // varible real key
    V(m,t+2*nr) =  1;
    
    // R
    V(m,t) =  r(t);
        
    //-X
    V(m,t+nr) = -x(t);

    // varible aimg key
    V(m+nb,t+2*nr+nr) = 1;
        
    //X
    V(m+nb,t) =  x(t);
    
    //R
    V(m+nb,t+nr) =  r(t);
    
    //find conected line
    up = find(destination==k);
    while (up <> [])

        //R
        V(m,up) =  r(up);
            
        //-X
        V(m,up+nr) = -x(up);
           
        //X
        V(m+nb,up) =  x(up);
        
        //R
        V(m+nb,up+nr) =  r(up);
        
        k =  source(up);
       
        //find conected line
        up = find(destination==k);
     end 
end

// Add Vel + A  
A = cat(1,A,V);


// --------------------------------------
// STEP 5 : Create Upper and Lower bounds
// ---------------------------------------
lei = 1;
// Upper Bound
ub = ones(nr*4,1) * vub;

// Lower Bounds
lb = ones(nr*4,1) * vlb;

for i = 1:nb 

    if (tipo(i) ~= 0) || (Pmax(i) ~= 0 || Pmin(i) ~=0) 

        pAr = zeros(nb*4,1);
        pAi = zeros(nb*4,1);

        pAr(i) =  1;
        pAi(i+nb) =  1;

        A = cat(2,A,pAr,pAi);       
        w = cat(1, w, 0.0, 0.0);
        gd = gd + 1;

        if (Qmax(i) ~= 0 || Qmin(i) ~=0)
            disp("Qmax and Qmin", "bus", [i],[Qmax(i)], [Qmin(i)]);
            ub(nr * 4 + 2*gd,1) = Qmax(i);
            lb(nr * 4 + 2*gd,1) = Qmin(i);
        end
        if (Pmax(i) ~= 0 || Pmin(i) ~=0)
            disp("Pmax and Pmin", "bus", [i],[Pmax(i)], [Pmin(i)]);
            ub((nr * 4 + 2*gd)-1,1) = Pmax(i);
            lb((nr * 4 + 2*gd)-1,1) = Pmin(i);
        end

        if (VM(i) ~= 0.0 && VM(i) ~= 1.00)            
            disp("Vr", "bus", [i],[VM(i)]);
            ub(nr * 2 + i,1) = VM(i);
            lb(nr * 2 + i,1) = VM(i);
        end    
    end
end


// --------------------------------------
// STEP 6 : Compute currente Load and update Ibus
// ---------------------------------------
//ZIP Model
for i=1 : nb
    
    if (Zr(i) ~= 0 || Zi(i) ~=0)
        Zc = Zr(i) + %i*Zi(i);
        Ibus(i) = Ibus(i) + Vf/Zc;
    end
    

    if (Ir(i) ~= 0 || Ii(i) ~=0)
        Ic = Ir(i) + %i*Ii(i);
        Ibus(i) = Ibus(i) + Ic;
    end
    

    if (Yr(i) ~= 0 || Yi(i) ~=0)
        Yc = Yr(i) + %i*Yi(i);
        Ibus(i) = Ibus(i) + Vf * Yc;
    end
 
 
    if (Pc(i) ~= 0 || Qc(i) ~=0)
        Sc = (Pc(i) + %i*Qc(i)) / Vf;
        Ibus(i) = Ibus(i) + Sc;
    end
    
 end

// b matrix = Ibus
b = cat(1, real(Ibus), imag(Ibus));  
b = cat(1,b,Vr_inicial);
b = cat(1,b,Vi_inicial);

//matriz H
H=eye(4*nr+(2*gd),4*nr+(2*gd));

// loss
sourceltaloss = 0;

//iteraction
iteraction = 0;

// sourcelta V inicial
sourceltaV = 1.2;

// Power P and Q 
S = Pesp + %i * Qesp;

Vm = 0;
       
tic(); 
// --------------------------------------
// STEP 9 : Convergence test while (Vp - V(p-1)) >= erro
// ---------------------------------------
while (sourceltaV >= erro)
    
    // last calculated Voltage
    lastV = abs(median(real(Vm)));

    // --------------------------------------   
    // STEP 7 : Solve the OCF problem
    // ---------------------------------------
    [I,iact,iter,f]=qpsolve(H,w,A,b,lb,ub,size(b)(1)-2)

    // currents km;
    Ikm =  I(1:nr) + %i * I(1+nr:nr+nr)
    
    // voltages;
    Vm = I(1+nr+nr:nr+nr+nr) + %i * I(1+nr+nr+nr:nr+nr+nr+nr);
    
    // curentes Ig
    Ig = I((4*nr)+1:4*nr+gd) - %i * I((4*nr)+1+gd:4*nr+1+gd);
    
    indices_barras = 1:length(Ikm);
    indices_gd =  1:length(Ig);
    
    disp("****************** Iteraction ********************"); 
    disp ([iteraction]);
    
    disp("****************** OCF - Currents km e Voltages Vk ********************");
    disp ([indices_barras', real(Ikm), imag(Ikm), real(Vm), imag(Vm)]);
    
    disp("****************** OCF - DERs current injections ********************"); 
    disp ([indices_gd', real(Ig), imag(Ig)]);
    
    disp("****************** OCF - sourcelta Vk ********************"); 
    disp ([sourceltaV]);

       
    // --------------------------------------
    // STEP 8 : Compute current Load and updatde Ibus
    // with the new Voltage
    // ---------------------------------------
    
    //new voltage
    Vm = cat(1,Vf,Vm);
   
    // new Ibus
    Ibus = zeros(nb,1);
    
    for i=1 : nb
        //Impedância constante
        if (Zr(i) ~= 0 || Zi(i) ~=0)
            Zc = Zr(i) + %i*Zi(i);
            Ibus(i) = Ibus(i) + Vm(i)/Zc;
        end
        
        //Corrente constante
        if (Ir(i) ~= 0 || Ii(i) ~=0)
            Ic = Ir(i) + %i*Ii(i);
            Ibus(i) = Ibus(i) + Ic;
        end
        
        //Admitância constante
        if (Yr(i) ~= 0 || Yi(i) ~=0)
            Yc = Yr(i) + %i*Yi(i);
            Ibus(i) = Ibus(i) + Vm(i)*Yc;
        end
 
        //Potência constante
        if (Pc(i) ~= 0 || Qc(i) ~=0)
            Sc = (Pc(i) + %i*Qc(i))/Vm(i);
            Ibus(i) = Ibus(i) + Sc;
        end
     end
              
    // refresh Ibus = b
    b = cat(1, real(Ibus), imag(Ibus));
    b = cat(1,b,Vr_inicial);
    b = cat(1,b,Vi_inicial);
    
    
    //sourceltaV  <= sourceltaVm(n) - sourceltaVm(n-1))
    sourceltaV = abs(median(real(Vm)) - lastV);
    
    //add iteraction counter
    iteraction = iteraction +1;
    
end

tempo_exec = toc();
disp(tempo_exec);

xlfont("SansSerif",10,%t,%t)
plot(real(Vm), "r")


title("Bus Voltage Profile 34-bus","fontsize",2);
xlabel('Bus');
ylabel('Vr[p.u.]');


