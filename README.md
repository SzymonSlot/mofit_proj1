# Projekt na zajęcia MOFiT "Dynamika punku materialnego"

Projekt zawiera implementację i analizę porównawczą różnych metod numerycznych rozwiązywania równań różniczkowych zwyczajnych (ODE) na przykładzie ruchu cząstki w złożonym polu potencjalnym.

## 🚀 Funkcje projektu

* **Implementacja wielu solverów**: Od prostego schematu Eulera i Verleta po zaawansowane metody RK4 i niejawne metody trapezów.
* **Adaptacyjny krok czasowy**: Algorytmy automatycznie dostosowujące gęstość próbkowania na podstawie zadanego poziomu tolerancji błędu.
* **Analiza fizyczna**: Monitorowanie zachowania energii całkowitej układu oraz generowanie portretów fazowych.

## 🔬 Model Fizyczny

Cząstka o masie $m$ porusza się w potencjale zdefiniowanym wzorem:
$$V(x) = -e^{-(x/l_1)^2} - 8e^{-(x-2)^2/l_2^2}$$
gdzie $l_1 = 1$ oraz $l_2 = 1/\sqrt{8}$. 

Układ uwzględnia również siłę oporu (tłumienie) zależną od prędkości: $F_{oporu} = -\alpha \cdot v$.

## 📁 Struktura plików

* `funkcje.py` – biblioteka zawierająca definicję potencjału, sił oraz główne algorytmy (`RK4`, `verlet`, `trapez`).
* `jawnyeuler.py`, `RK4.py`, `verlet.py` – skrypty do analizy zbieżności metod przy stałym kroku czasowym.
* `zadanie2.py` – porównanie metod adaptacyjnych (Euler vs Verlet vs RK4).
* `zadanie3.py` – analiza wpływu współczynnika tłumienia $\alpha$ na dobór kroku czasowego.
* `zadanie4.py` – implementacja niejawnej metody trapezów z wykorzystaniem metody Newtona.

## 🛠 Instalacja i uruchomienie

Wymagane biblioteki:
```
pip install numpy matplotlib numba pandas
```

📊 Przykładowe Wyniki

Program generuje zestaw czterech wykresów dla każdego przypadku:
1. Położenie $x(t)$ – trajektoria cząstki.
2. Prędkość $v(t)$ – zmiany prędkości w czasie.
3. Portret fazowy – zależność $v(x)$.
4. Energia $E(t)$ – weryfikacja stabilności numerycznej (zachowanie energii).
