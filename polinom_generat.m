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
    numar_total_variabile=na+nb; 
    combinatii=[]; 
    combinatie_curenta=zeros(1,numar_total_variabile); 

    function backtracking(combinatie_curenta, pozitie, suma_ramasa)

        if pozitie>numar_total_variabile
           
            if suma_ramasa>=0
                combinatii=[combinatii; combinatie_curenta];
            end
            return;
        end

        for i=0:suma_ramasa
            combinatie_curenta(pozitie)=i;
            backtracking(combinatie_curenta, pozitie+1, suma_ramasa-i);
        end
    end
  
    backtracking(combinatie_curenta, 1, m);
end