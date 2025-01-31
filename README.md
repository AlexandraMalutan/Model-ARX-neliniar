# Model-ARX-neliniar
Programarea unei functii care genereaza un model ARX neliniar de tip polinomial, pentru ordinele na, nb, si gradul m configurabile. De asemenea, programarea procdurii de regresie pentru identificarea parametrilor, si utilizarea modelului cu intrari noi.


Se da un set de date masurat de la un sistem dinamic cu o intrare si o iesire. Ordinul sistemului nu este mai mare de 3, si dinamica poate fi neliniara, ın timp ce iesirea poate fi afectata de zgomot. Vom dezvolta un model de tip cutie neagra pentru acest sistem, folosind o structura ARX neliniara de tip polinomial. Un al doilea set de date, masurat de la acelasi sistem, este furnizat pentru validarea modelului dezvoltat. Cele doua seturi de date sunt furnizate ıntr-un fisier MATLAB, respectiv ın variabilele id si val, ambele obiecte de tip iddata din toolbox-ul de identificarea sistemelor. Reamintim ca intrarea, iesirea, si perioada de esantionare sunt disponibile in campurile u, y, Ts ale acestor obiecte. In caz ca toolbox-ul nu este instalat, ˘ acelesi seturi de date sunt furnizate si in format vectorial, id array si val array, fiecare dintre ele o matrice cu structura: valorile de timp pe prima coloana, intrarea pe a doua, si iesirea pe ultima coloana.
