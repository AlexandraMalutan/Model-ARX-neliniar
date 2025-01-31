% incarcarea fisierelor cu datele primite
load('iddata-11.mat')
figure;
plot(id.InputData{1}), title("Intrarea sistemului pentru datele de identificare")
figure;
plot(id.OutputData{1}), title("Iesirea sistemului pentru datele de identificare")

figure;
plot(val.InputData{1}), title("Intrarea sistemului pentru datele de validare")
figure;
plot(val.OutputData{1}), title("Iesirea sistemului pentru datele de validare")

na=2;
nb=2;
nk=1;
m=4;

idU=id.InputData{1};
idu=detrend(idU);
valU=val.InputData{1};
valu=detrend(valU);
idY=id.OutputData{1};
idy=detrend(idY);
valY=val.OutputData{1};
valy=detrend(valY);


valoaremseid=zeros(1,m);
valorimseval=zeros(1,m);
valoaremsesimidentificare=zeros(1,m);
valoaremsesimvalidare=zeros(1,m);
%pentru simulare
Nid=length(idy); %dimensiunea setului de date de identificare
ysimulat_id=zeros(Nid, 1); %initializez vectorul pentru iesirile simulate, cu valori 0
ysimulat_id(1:na)=idy(1:na); %realizarea conditiilor initiale prin setarea valorilor simulate la valorile reale
Nval=length(valy);%dimensiunea setului de date de validare
ysimulat_val=zeros(Nval,1);%initializez vectorul pentru iesirile simulate, cu valori 0
ysimulat_val(1:na)=valy(1:na);



for l=1:m
    %pentru identificare
    phi=polinom_generat2(idu,idy,l,na,nb,nk);%Construieste matricea de regresie pentru identificare
    theta=phi\idy;% Calculeaza parametrii modelului folosind regresia liniara
    yaproximat=phi*theta; %calculeaza iesirile aproximate pe datele de identificare
    %Calculeaza eroarea medie patratica pe setul de identificare
    suma=0;
    for i=1:length(idu)
        suma=suma+(yaproximat(i)-idy(i)).^2;
    end
    mseid=(1/length(idu))*suma;
    valoaremseid(l)=mseid;

    %pentru validare
    phivalidare=polinom_generat2(valu,valy,l,na,nb,nk);% Construieste matricea de regresie pentru validare
    yaproximatvalidare=phivalidare*theta;%Calculez iesirile aproximate pentru datele de validare
    %Calculez eroarea medie patratica pe setul de validare
    suma=0;
    for i=1:length(valu)
        suma=suma+(yaproximatvalidare(i)-valy(i)).^2;
    end
    mseval=(1/length(valu))*suma;
    valorimseval(l)=mseval;

    phi_simulat=zeros(Nid, size(phi, 2));
    %pentru simularea identificarii
    for i=1+na:Nid  %aici se itereaza de la valoarea na pentru fiecare pas 
        phi_simulat=polinom_generat2(idu(1:i), ysimulat_id(1:i), l, na, nb, nk);%aici noi construim phi pentru pasul la care ne aflam
        ysimulat_id(i)=phi_simulat(end, :)*theta; % i se atribuie iesirii simualte valoarea calculata anterior din phi pentru pasul 
        % curent inmultit cu vectorul theta pentru realizarea simularii
    end
    suma=0;
    for i=1:Nid
        suma=suma+(ysimulat_id(i)-idy(i)).^2;
    end
    msesimidentificare=(1/Nid)*suma;
    valoaremsesimidentificare(l)=msesimidentificare;

    phi_simulat_val=zeros(Nval, size(phivalidare, 2));%initioalizez si matricea de creare a matricei phi pentru simularea validarii
    for i=1+na:Nval
    phi_simulat_val=polinom_generat2(valu(1:i), ysimulat_val(1:i), l, na, nb, nk);
    ysimulat_val(i)=phi_simulat_val(end,:)*theta;
    end
    suma=0;
    for i=1:Nval
        suma=suma+(ysimulat_val(i)-valy(i)).^2;
    end
    msesimvalidare=(1/Nval)*suma;
    valoaremsesimvalidare(l)=msesimvalidare;
end


disp("Eroare MSE pe datele de identificare:")
disp(mseid)
figure
plot(1:m, valoaremseid), title("MSE pentru identificare"); 

disp("Eroare MSE pe datele de validare:")
disp(mseval)
figure;
plot(1:m,valorimseval), title("MSE pentru validare");

disp("Eroare MSE pentru simularea pe datele de identificare:")
disp(msesimidentificare)
figure;
plot(1:m,valoaremsesimidentificare), title("MSE pe simulare identificare");


disp("Eroare MSE pentru simularea pe datele de validare:")
disp(msesimvalidare)
figure;
plot(1:m,valoaremsesimvalidare), title("MSE pe simulare validare");

figure;
plot(idy)
hold on;
plot(yaproximat, 'r-')
legend('Iesirea reala','Iesirea modelata')
title('Graficul pentru identificare');

figure;
plot(valy)
hold on;
plot(yaproximatvalidare, 'r-')
legend('Iesirea sistemului','Iesirea modelului')
title('Graficul pentru validare');

figure;
plot(idy)
hold on;
plot(ysimulat_id, 'r-')
legend('Iesirea modelului','Iesirea simulata')
title('Simulare pe identificare');

figure;
plot(valy)
hold on;
plot(ysimulat_val, 'r-')
legend('Iesirea modelului','Iesirea simulata')
title('simulare pentru validare');


























function phi = polinom_generat2(x, y, m, na, nb, nk)
    dim=length(y); %dimensiunea setului de date
    combinatii=generare_combinatii(na, nb, m);%generarea tuturor combinatiilor posibile
    %de exponenti pentru polinom ce il creem
    nrtermeni=size(combinatii, 1);%nr de termeni din polinom, reprezintand elementele 
    %de pe coloana din matricea phi
    phi=zeros(dim, nrtermeni);

    %incepem sa construim matricea phi, pe fiecare rand
    for i=1:dim
        %parcurgem toate combinatiile de exponenti, stiind ca puterea nu
        %trebuie sa depaseasca gradul m al polinomului 
        for j=1:nrtermeni
            exponent=combinatii(j, :);%expoentii pe linia curenta pentru toate coloanele
            termen=1;%initial termenul este 1

            %calculam termenii pentru iesirile intarziate pana la ordinul na
            for k=1:na
                if (i-k>0)
                    termen=termen*(-y(i-k)^exponent(k));
                else
                    termen=0;
                end
            end

            %calculez termenii pentru intrarile intarziate
            for k=1:nb
                if (i-nk-k+1>0)
                   termen=termen*(x(i-nk-k+1)^exponent(na+k));
                else
                   termen=0;
                end
            end
            %temrenul calculat se introduce in patricea phi la pozitia
            %corespunzoare
            phi(i, j)=termen;
         end
     end
end






function combinatii=generare_combinatii(na, nb, m)
    numar_total_variabile=na+nb; %suma iesirilor intarziate si a intrarilor intarziate
    combinatii=[]; %combinatiile valide sunt stocate in acest vector 
    combinatie_curenta=zeros(1,numar_total_variabile); %acesta este un vector ajutator fiind o combinatie temporare
    %aceasta se schimba la fiecare pas (combinatia curent de exponenti)

    %Am folosit functia recursiva prin metoda backtracking pentru generarea
    %combinatiilor
    %pozitie este pozitia curenta in vectorul de exponenti
    %suma_ramasa este gradul maxim ramas de alocat intre variabile
    function backtracking(combinatie_curenta, pozitie, suma_ramasa)

        if pozitie>numar_total_variabile
            %daca gradul maxim ramas este 0 sau mai mare salvam combinatia 
            if suma_ramasa>=0
                combinatii=[combinatii; combinatie_curenta]; %adaug combinatia valida
            end
            return;%apelul recursiv care indica iesirea din functie
        end

        %parcurg toate valorile posibile pentru exponentul variabilei
        %curente
        for i=0:suma_ramasa
            combinatie_curenta(pozitie)=i; %setez exponentul 
            %apelez recursiv functia pentru urmatoarea variabila
            backtracking(combinatie_curenta, pozitie+1, suma_ramasa-i);
        end
    end
    %realizez acest proces de la prima variabila pana cand se parcurg toate
    %combinatiile si ajung la polinom rezultat pentru gradul dat m
    backtracking(combinatie_curenta, 1, m);
end