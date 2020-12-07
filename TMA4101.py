'''
Denne modulen er laget av Sebastian W. Kalland til 
TMA4101 høsten 2020, ideelt sett for lettere gjennomføring av eksamen,
men kan like så godt brukes til senere formål også.

Den burde ligge under: C:\ProgramData\Anaconda3\Lib\TMA4101.py
da er den tilgjengelig for import i Spyder overalt på maskinen.

Versjon 03.12.2020
'''

def newtons_metode(f, x_0, n, **kwargs):
    '''
    Newtons metode
    ----------
    Tar inn en initialgjetning, et funksjonsuttrykk, og dens deriverte.
    Gir deretter et nullpunkt ut i fra initialgjetningen.
    
    
    Parameters
    ----------
    f : Funksjon
        Den matematiske funksjonen vi skal finne nullpunkt til
    x_0 : Float
        Initialgjetningen vi starter metoden på
    n : Integer
        Antall kjøringer før resultatet skal returneres
    
    
    Kwargs
    ----------
    p : True/False
        Printer per iterasjon om satt til True
    df : Funksjon
        Den matematiske funksjonens deriverte dersom tilgjengelig, 
        uten verdi vil en approksimasjon i punktet gjennomføres
    form : String
        Formatet printen skal komme ut i. Default: "15.10f"
    h : Float
        Dersom det må gjøres en numerisk utregning på den deriverte er dette
        stegstørrelsen for gjetningen. Se numerisk_derivasjon() 
        for mer informasjon. Default: 0.001


    Returns
    -------
    x-verdi etter n iterasjoner
    '''
    x = x_0
    
    h = kwargs.get("h") if kwargs.get("h") else 0.001
    
    form = kwargs.get("form") if kwargs.get("form") else "15.10f"
    
    for i in range(n):
        if kwargs.get("df")==None:
            df = numerisk_derivasjon(f,x,h=h)
            x = x - f(x) / df
        else:
            x = x - f(x) / kwargs.get("df")(x)
            
        if kwargs.get("p"):
            print("i:" + format(i,"5d") + "\t\t x: " + format(x,form))
    
    if kwargs.get("p"):
        print("Sluttresultatet ble: x = " + format(x,form))
    
    return x

def fikspunktiterasjon(f, x_0, n, **kwargs):
    '''
    Fikspunktiterasjon
    ----------
    Tar inn en initialgjetning, og et funksjonsuttrykk som skal være lik x 
    Fungerer kun når |f'(x)| < 1
    
    
    Parameters
    ----------
    f : Funksjon
        Den matematiske funksjonen vi skal finne som skal bli lik x
    x_0 : Float
        Initialgjetningen vi starter iterasjonen med
    n : Integer
        Antall iterasjoner før resultatet skal returneres
    
    
    Kwargs
    ----------
    p : True/False
        Printer per iterasjon om satt til True
    form : String
        Formatet printen skal komme ut i. Default: "15.10f"


    Returns
    -------
    x-verdi etter n iterasjoner
    '''
    x = x_0
    
    form = kwargs.get("form") if kwargs.get("form") else "15.10f"
    
    for i in range(n):
        x = f(x)
            
        if kwargs.get("p"):
            print("i:" + format(i,"5d") + "\t\t x: " + format(x,form))
    
    if kwargs.get("p"):
        print("Sluttresultatet ble: x = " + format(x,form))
    
    return x

def numerisk_derivasjon(f, x, **kwargs):
    '''
    Numerisk derivasjon
    ----------
    Tar inn en funksjon og en x-verdie funksjonens deriverte skal regnes ut i.
    Bruker formelen som har taylorrekkens ledd med funksjonens fjerdederiverte som
    feilmargin:
        df = (f(x - 2*h) - 8*f(x - h) + 8*f(x + h) - f(x + 2*h)) / (12*h)
    
    
    Parameters
    ----------
    f : Funksjon
        Funksjonen som skal numerisk deriveres i punktet x
    x : Float
        Punktet på x-aksen funksjonen skal deriveres i
        
    
    Kwargs
    ----------
    h : Float
        Hvor store steg som skal tas i funksjonens utregning av den deriverte, 
        generelt sett er det slik at jo mindre jo bedre. Default: 0.001


    Returns
    -------
    df : Float
        Den deriverte av f i punktet x, f(x)
    '''
    h = kwargs.get("h") if kwargs.get("h") else 0.001
    
    df = f(x - 2*h) - 8*f(x - h) + 8*f(x + h) - f(x + 2*h)
    df /= 12*h
    
    return df

def deriver(f, x):
    '''
    Deriver
    ----------
    Returnerer funksjonen f sin deriverte i f, f(x).
    Se numerisk_derivasjon.
    
    
    Parameters
    ----------
    f : Funksjon
        Funksjonen f
    x : Float
        x-verdien den deriverte skal hentes i
        
    
    Returns
    ----------
    Returnerer "numerisk_derivasjon" med h = 0.001
    '''
    return numerisk_derivasjon(f, x)

def venstre_riemannsum(f, a, b, n):
    '''
    Venstre_riemannsum
    ----------
    Returnerer funksjonen f sin riemannsum på intervallet [a,b]
    

    Parameters
    ----------
    f : Funksjon
        Den matematiske funksjonen vi skal finne riemannsummen under
    a : Float
        Startsverdi av riemannsummen. Den verdien lengst til venstre.
    b : Float
        Sluttverdi til riemannsummen. Den verdien lengst til høyre.
    n : Int
        Antall punkter intervallet [a,b] skal deles opp i og vurderes for.


    Returns
    -------
    Float
        Verdien av integralet for funksjonen f på intervallet [a,b] regnet ut ved
        venstre riemannsum med n antall partisjoner.
    '''
    s = 0
    
    for i in range(n):
        x = a + (b - a) * i/n
        s += f(x) * (abs(a - b)/n)
    
    return s

def høyre_riemannsum(f, a, b, n):
    '''
    Høyre_riemannsum
    ---------
    Se venstre_riemannsum.
    
    
    Returns
    ---------
    Float
        venstre_riemannsum evaluert ett punkt mer til høyre 
        for alle punkter i partisjonen
    '''
    return venstre_riemannsum(f, 
                              a + abs(a-b)/n, 
                              b + abs(a-b)/n, 
                              n)

def riemannsum(f, a, b, n):
    '''
    Riemannsum
    ---------
    Se venstre_riemannsum.
    
    
    Returns
    ---------
    Float
        venstre_riemannsum evaluert for en partisjon på n antall punkter.
    '''
    return venstre_riemannsum(f, a, b, n)

def trapesmetoden_riemann_vAvg(f, a, b, n):
    '''
    Trapesmetoden_riemann_vAvg
    ---------
    Tar gjennomsnittet av venstre- og høyre riemannsum.
    Se venstre_riemannsum og høyre_riemannsum.
    
    
    Returns
    ---------
    Float
        venstre_riemannsum og høyre_riemannsum sitt gjennomsnitt.
    '''
    return venstre_riemannsum(f, a, b, n)/2 + høyre_riemannsum(f, a, b, n)/2

def trapesmetoden_riemann(f, a, b, n):
    '''
    Trapesmetoden_riemann
    ---------
    f(a)/(2*n) + venstre_riemannsum uten første ledd + f(b)/(2*n)
    Se venstre_riemannsum.
    
    
    Returns
    ---------
    Float
    '''
    return f(a)/(2*n) + venstre_riemannsum(f, a+(abs(a-b)/n), b, n-1) + f(b)/(2*n)

def midtpunktregelen(f, a, b, n):
    '''
    Midtpunktregelen
    ---------
    Tar riemannsum evaluert i punktene midt mellom punktene venstre-
    og høyre riemannsum evalueres i.
    Se venstre_riemannsum.
    
    
    Returns
    ---------
    Float
        en riemannsum evaluert for 1/2 punkt lengre til høyre enn venstre riemannsum.
    '''
    return venstre_riemannsum(f, 
                              a + abs(a-b)/(2*n), 
                              b + abs(a-b)/(2*n), 
                              n)
    
def integralet(f, a, b):
    '''
    Integralet
    ---------
    Se trapesmetoden_riemann.
    
    
    Returns
    ---------
    Float
        trapesmetoden på f på [a,b] evaluert for en partisjon på 100 000 punkter.
    '''
    return trapesmetoden_riemann(f, a, b, 100_000)

def euler_eksplisitt(dy, x_0, y_0, n, x_n, **kwargs):
    '''
    Euler eksplisitt
    ----------
    Tar et uttrykk for den deriverte av y og x og en startsverdi, og lager en graf etter de
    betingelsene. Går som regel "utenfor" grafens bue, hvorav implisitt går "innenfor".
    Merk at dy må ta inn x og y i dy, selv om ikke begge brukes.
    
    
    Parameters
    ----------
    dy : Funksjon
        Funksjonsuttrykk for y' av x OG y.
        Eksempelvis y' = y^2 - 3x, da er dy = lambda x, y: y**2 - 3*x
    x_0 : Float
        Første x-verdi
    y_0 : Float
        Første y-verdi
    n : Int
        Antall oppdelinger
    x_n : Float
        Siste x-verdi
    
    
    Kwargs
    ---------
    plot : True/False
        Om funksjonen skal plotte en graf. Default True.
    implisitt : True/False
        Å heller gjøre euler implisitt. Default False.
    implisitt_n : Int
        Hvor nøyaktig den implisitte gjetningen er. Default 64.
    show : True/False
        Plot som egen graf, eller tilhørende andre grafer? 
        Default False, ergo tilhørende andre grafer.

    Returns
    -------
    [x, y]: Array
        numpy arrays for x og tilhørende y i et eget numpy array
    '''
    from numpy import zeros, linspace, array
    from matplotlib.pyplot import plot, xlabel, ylabel, show, legend
    
    impl = True if kwargs.get("implisitt") else False
    impl_n = kwargs.get("implisitt_n") if kwargs.get("implisitt_n") else 64
    
    h = abs(x_0-x_n)/n #Steglengde
    
    y = zeros(n)
    y[0] = y_0
    
    x = linspace(x_0, x_n, num=n)
    
    for i in range(n-1):
        y[i+1] = y[i] + h*dy(x[i], y[i])
        
        if impl:
            for _ in range(impl_n): #Approksimerer neste y ved å se på neste y's deriverte og justerer mange ganger så det ser rett ut
                y[i+1] = y[i] + h*dy(x[i + 1], y[i + 1])
    
    if kwargs.get("plot")!=False: #Plotter by default
        plot(x, y, label=("Euler " + ("implisitt" if impl else "eksplisitt")))
        xlabel("x")
        ylabel("y")
        legend()
        if kwargs.get("show"): show()
    
    return array([x, y])

def euler_implisitt(dy, x_0, y_0, n, x_n, **kwargs):
    '''
    Parameters
    ----------
    Se euler_eksplisitt

    Returns
    -------
    euler_eksplisitt(dy, x_0, y_0, n, x_n, implisitt=True, **kwargs)

    '''
    try:
        del kwargs["implisitt"]
    except:
        pass
    return euler_eksplisitt(dy, x_0, y_0, n, x_n, implisitt=True, **kwargs)

def trapesmetoden_difflikning(dy, x_0, y_0, n, x_n, **kwargs):
    '''
    Trapesmetoden_difflikning
    ----------
    Tar et uttrykk for den deriverte av y og x og en startsverdi, og lager en graf etter de
    betingelsene. Tar gjennomsnittsverdien av euler eksplisitt sitt steg og implisitt sitt steg
    som sin neste y-verdi.
    Merk at dy må ta inn x og y i dy, selv om ikke begge brukes.
    
    
    Parameters
    ----------
    dy : Funksjon
        Funksjonsuttrykk for y' av x OG y.
        Eksempelvis y' = y^2 - 3x, da er dy = lambda x, y: y**2 - 3*x
    x_0 : Float
        Første x-verdi
    y_0 : Float
        Første y-verdi
    n : Int
        Antall oppdelinger
    x_n : Float
        Siste x-verdi
    
    
    Kwargs
    ---------
    plot : True/False
        Om funksjonen skal plotte en graf. Default True.
    implisitt_n : Int
        Hvor nøyaktig den implisitte gjetningen er. Default 64.
    show : True/False
        Plot som egen graf, eller tilhørende andre grafer?

    Returns
    -------
    x, y: Array
        numpy arrays for x og tilhørende y.
    '''
    from numpy import zeros, linspace
    from matplotlib.pyplot import plot, xlabel, ylabel, show, legend
    
    impl_n = kwargs.get("implisitt_n") if kwargs.get("implisitt_n") else 64
    
    h = abs(x_0-x_n)/n #Steglengde
    
    y = zeros(n)
    y[0] = y_0
    
    x = linspace(x_0, x_n, num=n)
    
    for i in range(n-1):
        y[i+1] = y[i] + h*dy(x[i], y[i])
        for _ in range(impl_n): #Approksimerer neste y ved å se på neste y's deriverte og justerer mange ganger så det ser rett ut
            y[i+1] = y[i] + h*dy(x[i + 1], y[i + 1])
        y[i+1] = y[i] + h * (dy(x[i], y[i]) + dy(x[i+1], y[i+1]))/2
    
    if kwargs.get("plot")!=False: #Plotter by default
        plot(x, y, label="Trapesmetoden")
        xlabel("x")
        ylabel("y")
        legend()
        if kwargs.get("show"): show()
    
    return x, y

def heuns_metode(dy, x_0, y_0, n, x_n, **kwargs):
    '''
    Heuns_metode
    ----------
    Tar et uttrykk for den deriverte av y og x og en startsverdi, og lager en graf etter de
    betingelsene. Tar gjennomsnittsverdien av euler eksplisitt sitt steg og implisitt steg
    vurdert i eksplisitt y.
    Merk at dy må ta inn x og y, selv om ikke begge brukes.
    
    
    Parameters
    ----------
    dy : Funksjon
        Funksjonsuttrykk for y' av x OG y.
        Eksempelvis y' = y^2 - 3x, da er dy = lambda x, y: y**2 - 3*x
    x_0 : Float
        Første x-verdi
    y_0 : Float
        Første y-verdi
    n : Int
        Antall oppdelinger
    x_n : Float
        Siste x-verdi
    
    
    Kwargs
    ---------
    plot : True/False
        Om funksjonen skal plotte en graf. Default True.
    show : True/False
        Plot som egen graf, eller tilhørende andre grafer?

    Returns
    -------
    x, y: Array
        numpy arrays for x og tilhørende y.
    '''
    from numpy import zeros, linspace
    from matplotlib.pyplot import plot, xlabel, ylabel, show, legend
    
    h = abs(x_0-x_n)/n #Steglengde
    
    y = zeros(n)
    y[0] = y_0
    
    x = linspace(x_0, x_n, num=n)
    
    for i in range(n-1):
        y[i+1] = y[i] + h*dy(x[i], y[i])
        y[i+1] = y[i] + h * (dy(x[i], y[i]) + dy(x[i+1], y[i+1]))/2
    
    if kwargs.get("plot")!=False: #Plotter by default
        plot(x, y, label="Heuns metode")
        xlabel("x")
        ylabel("y")
        legend()
        if kwargs.get("show"): show()
    
    return x, y

def midtpunktmetoden(dy, x_0, y_0, n, x_n, **kwargs):
    '''
    Midtpunktmetoden
    ----------
    Tar et uttrykk for den deriverte av y og x og en startsverdi, og lager en graf etter de
    betingelsene. Finner neste y av euler implisitt, og setter så neste y til å være evaluert for
    midten av dette stegets x- og y-verdier og neste stegs verdier.
    Merk at dy må ta inn x og y i dy, selv om ikke begge brukes.
    
    
    Parameters
    ----------
    dy : Funksjon
        Funksjonsuttrykk for y' av x OG y.
        Eksempelvis y' = y^2 - 3x, da er dy = lambda x, y: y**2 - 3*x
    x_0 : Float
        Første x-verdi
    y_0 : Float
        Første y-verdi
    n : Int
        Antall oppdelinger
    x_n : Float
        Siste x-verdi
    
    
    Kwargs
    ---------
    plot : True/False
        Om funksjonen skal plotte en graf. Default True.
    implisitt_n : Int
        Hvor nøyaktig den implisitte gjetningen er. Default 64.
    show : True/False
        Plot som egen graf, eller tilhørende andre grafer?

    Returns
    -------
    x, y: Array
        numpy arrays for x og tilhørende y.
    '''
    from numpy import zeros, linspace
    from matplotlib.pyplot import plot, xlabel, ylabel, show, legend
    
    impl_n = kwargs.get("implisitt_n") if kwargs.get("implisitt_n") else 64
    
    h = abs(x_0-x_n)/n #Steglengde
    
    y = zeros(n)
    y[0] = y_0
    
    x = linspace(x_0, x_n, num=n)
    
    for i in range(n-1):
        y[i+1] = y[i] + h*dy(x[i], y[i])
        for _ in range(impl_n): #Approksimerer neste y ved å se på neste y's deriverte og justerer mange ganger så det ser rett ut
            y[i+1] = y[i] + h*dy(x[i + 1], y[i + 1])
        y[i+1] = y[i] + h * dy((x[i]+x[i+1])/2, (y[i]+y[i+1])/2)
    
    if kwargs.get("plot")!=False: #Plotter by default
        plot(x, y, label="Trapesmetoden")
        xlabel("x")
        ylabel("y")
        legend()
        if kwargs.get("show"): show()
    
    return x, y

def matrisemultiplikasjon(a,b):
    '''
    Bruker numpy til å regne ut produktet av matrisene

    Parameters
    ----------
    a : List/Array/Matrix
        Venstre matrise
    b : List/Array/Matrix
        Høyre matrise

    Returns
    -------
    Array
        Matrisenes produkt

    '''
    from numpy import array, matmul
    a = array(a)
    b = array(b)
    return matmul(a,b)

def prikkprodukt(a,b):
    '''
    Bruker numpy til å regne ut prikkproduktet av to vektorer

    Parameters
    ----------
    a : List/Array
        DESCRIPTION.
    b : List/Array
        DESCRIPTION.

    Returns
    -------
    Array
        prikkproduktet av arrayene

    '''
    from numpy import array, dot
    a = array(a)
    b = array(b)
    return dot(a,b)
    
def determinant(a):
    '''
    Bruker numpy.linalg til å regne ut determinanten til matrisen a

    Parameters
    ----------
    a : Array/Matrix

    Returns
    -------
    Float
        Determinanten til a

    '''
    from numpy.linalg import det
    return det(a)

def kryssprodukt(a,b):
    '''
    Bruker numpy til å regne ut kryssproduktet til vektorene a og b

    Parameters
    ----------
    a : List/Array
        Venstre vektor
    b : List/Array
        Høyre vektor

    Returns
    -------
    Float
        Kryssproduktet a x b

    '''
    from numpy import array, cross
    a = array(a)
    b = array(b)
    return cross(a,b)

def buelengde(df,a,b,n):
    '''
    Bruker trapesmetoden til å numerisk regne ut buelengden av en funksjon

    Parameters
    ----------
    df : Funksjon
        Den deriverte av funksjonen du ønsker buelengden til
    a : Flaot
        Startverdi for x
    b : Float
        Sluttverdi for x
    n : Int
        Oppløsning til trapesmetoden

    Returns
    -------
    Float
        Buelengden fra a til b på f

    '''
    from numpy import sqrt
    f = lambda x: sqrt(1 + df(x)**2)
    return trapesmetoden_riemann(f, a, b, n)

def sylinderskall(f,a,b,n):
    '''
    Bruker trapesmetoden til å numerisk regne ut volumet av et sylinderskall med høyde f fra radius a til b

    Parameters
    ----------
    f : Funksjon
        Funksjonen for høyden
    a : Float
        Indre radius
    b : Float
        Ytre radius
    n : Int
        Oppløsning til trapesmetoden

    Returns
    -------
    Float
        Volumet av sylinderskallet

    '''
    from numpy import pi
    g = lambda x: x * f(x)
    return 2 * pi * trapesmetoden_riemann(g, a, b, n)

def dreievolum(f,a,b,n):
    '''
    Bruker trapesmetoden til å numerisk regne ut volumet av en funksjon for radius rundt x-aksen fra a til b

    Parameters
    ----------
    f : Funksjon
        Funksjon for radius
    a : Float
        Startsverdi for x
    b : Float
        Sluttverdi for x
    n : Int
        Oppløsning til trapesmetoden

    Returns
    -------
    Float
        Volumet av dreieelementet

    '''
    from numpy import pi
    g = lambda x: f(x)*f(x)
    return pi * trapesmetoden_riemann(g, a, b, n)

if __name__ == "__main__":
    print('''------
Det skal nå kjøres en rekke tester for å vise at modulen fungerer.    
------''')
    f = lambda x: x**2 - x
    print("newtons metode av f(x) = x^2 - x fra 3 med 7 iterasjoner burde bli 1, og ble:", newtons_metode(f, 3, 7))
    print("numerisk derivasjon av f(x) = x^2 - x i x = 3 burde bli 5, og ble:", deriver(f,3))
    f = lambda x: x/2 + 1/2
    print("fikspunktiterasjon av f(x) = x/2 + 1/2 fra 1 med 6 iterasjoner burde bli 1, og ble:", fikspunktiterasjon(f, 1, 6))
    f = lambda x: x
    print("trapesmetoden (integral) for integralet under f(x) = x fra 0 til 2 burde gi 2, og gav:", trapesmetoden_riemann_vAvg(f, 0, 2, 1000))
    dy = lambda x,y: x**2 + y**2
    trapesmetoden_difflikning(dy, 0, 0, 1000, 1)
    print("Sjekk plots for å se at Trapesmetodens (difflikning) graf ble plottet for y'=x^2+y^2")
    print("[[1,2],[3,4]] ganger [[5,6],[7,8]] burde bli [[19,22],[43,50]], det ble:",matrisemultiplikasjon([[1,2],[3,4]], [[5,6],[7,8]]))
    f = lambda x: x
    print("Dreievolumet for f(x)=x fra 0 til 1 burde være 1/3 * pi, og det er:", dreievolum(f, 0, 1, 1000))
    f = lambda x: -x + 1
    print("Sylinderskallet for f(x)=-x + 1 fra 0 til 1 burde være 1/3 pi, og det er", sylinderskall(f, 0, 1, 1000))
    df = lambda x: 2*x
    print("Buelengden for f(x)=x^2 fra 0 til 1 burde være lengre enn 0.73, og kortere enn 2:", buelengde(df, 0, 1, 1000))