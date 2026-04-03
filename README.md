# Metoda Elementów Skończonych – sprawozdanie końcowe

**Grupa lab. 6**  
**Data:** 14.01.2026  
**Imię i nazwisko:** Łukasz Żegliński  
**Temat:** Metoda Elementów Skończonych – sprawozdanie końcowe

## 1. Wstęp i informacje o rozwiązywanym problemie

### a) Temat sprawozdania

Tematem sprawozdania jest opis implementacji oraz działania programu, mającego na celu symulację niestacjonarnej wymiany ciepła. Program dokonuje obliczeń przy wykorzystaniu Metody Elementów Skończonych.

### b) Symulacja problemu rozkładu temperatury – opis teoretyczny

Celem programu jest wyznaczenie temperatur w węzłach elementów na zdefiniowanym obszarze. Odbywa się to poprzez rozwiązanie równania różniczkowego w kontekście niestacjonarnej wymiany ciepła. W tym celu zastosowano III warunek brzegowy (konwekcyjny):

\[q = \alpha(t - t_{ot})\]

gdzie:  
1. \(q\) – gęstość strumienia ciepła  
2. \(\alpha\) – współczynnik przejmowania ciepła  
3. \(t\) – temperatura na powierzchni ciała  
4. \(t_{ot}\) – temperatura otoczenia

Warunek brzegowy jest reprezentowany dzięki zastosowaniu macierzy \(H_{bc}\):

\[[H_{BC}] = \int_{S} \alpha(\{N\}\{N\}^{T}) dS = \sum_{i=1}^{n_{pc}} f(pc_i) w_i \det[J]\]

gdzie:  
- \(H_{bc}\) – macierz przedstawiająca warunek brzegowy dla danego elementu  
- \(\int_S\) – całka po powierzchni elementu  
- \(\alpha\) – współczynnik przejmowania ciepła  
- \(N\) – wektor funkcji kształtu dla danego elementu  
- \(\sum_{i=1}^{n_{pc}}\) – suma po punktach całkowania  
- \(w_i\) – waga zależna od schematu całkowania  
- \(\det[J]\) – wyznacznik jakobianu przekształcenia

\[[P] = \int_S \alpha\{N\} t_{ot} dS = \sum_{i=1}^{n_{pc}} f(pc_i) w_i \det[J]\]

gdzie:  
- \(P\) – wektor obciążeń  
- \(f\) – całka po powierzchni elementu  
- \(\alpha\) – współczynnik przejmowania ciepła  
- \(N\) – wektor funkcji kształtu dla danego elementu  
- \(t_{ot}\) – temperatura otoczenia  
- \(w_i\) – waga zależna od schematu całkowania  
- \(\det[J]\) – wyznacznik jakobianu przekształcenia


## 2. Charakterystyka oprogramowania

Program został zaimplementowany w języku C++.

### a) Działanie programu krok po kroku:

1. Wczytanie siatki do obiektu `Grid`.
2. Utworzenie obiektu typu `GlobalMatrix` – jego zadaniem jest przechowywanie globalnych macierzy \(C\), \(H+bc\) oraz wektora \(P\).
3. Dla każdego punktu w siatce:  
   I. Utworzenie obiektu `ElementData`  
   II. Obliczenie lokalnej macierzy \(H\)  
   III. Obliczenie lokalnej macierzy \(bc\) i wektora \(P\) (operacje wykonywane są na raz w jednej funkcji)  
   IV. Dodanie lokalnych macierzy \(H\) i \(bc\)  
   V. Obliczenie lokalnej macierzy \(C\)  
   VI. Dodanie lokalnych danych do struktury `GlobalMatrix`
4. Dla każdego kroku czasowego rozwiązanie układu równań z zastosowaniem rozkładu na macierze \(L\) i \(U\) (algorytm Doolittle’a).
5. Zapis wyników do pliku.

### b) Odpowiedzialności klas:

- **Node** – przechowywanie położenia punktu całkowania i informacji o warunku brzegowym.
- **Element** – przechowuje ID przypisanych node'ów.
- **Jakobian** – przechowywanie i obliczanie jakobianu przekształcenia.
- **Grid** – przechowywanie informacji o siatce (tablica node, tablica element, dane materiałowe: Density, conductivity, SpecificHeat itp.).
- **ElementData** – klasa odpowiedzialna za obliczanie lokalnych macierzy dla node'a.
- **GlobalMatrix** – przechowuje i uzupełnia globalne macierze \(C\), \(H_{bc}\) i wektor \(P\).

## 3. Porównanie wyników i wizualizacja

Testy oprogramowania przeprowadzono na podstawie trzech siatek dostarczonych przez prowadzącego (`test1.txt`, `MixGrid.txt`, `31_31.txt`) przy użyciu 2, 3 i 4-punktowego schematu całkowania. Kryterium porównawczym są wyniki min i max temperatur w węzłach dla każdej iteracji, przedstawione przez prowadzącego.

### a) Siatka 4x4
| Schemat 3 pkt: min | max | Różnica min [%] | Różnica max [%] | Schemat 4 pkt: min | max | Różnica min [%] | Różnica max [%] |
|--------------------|-----|-----------------|-----------------|--------------------|-----|-----------------|-----------------|
| 110,038 | 365,815 | 0,00000021 | -0,00000129 | 110,038 | 365,815 | 0,00000021 | -1,00000129 |
| 168,837 | 502,592 | -0,00000010 | 0,00000057 | 168,837 | 502,592 | -0,00000010 | -1,00000057 |
| 242,801 | 587,373 | 0,00000060 | 0,00000057 | 242,801 | 587,373 | 0,00000060 | -1,00000057 |
| 318,615 | 649,387 | 0,00000128 | -0,00000074 | 318,615 | 649,387 | 0,00000128 | -2,00000074 |
| 391,256 | 700,068 | 0,00000053 | -0,00000060 | 391,256 | 700,068 | 0,00000053 | -3,00000060 |
| 459,037 | 744,063 | 0,00000021 | -0,00000046 | 459,037 | 744,063 | 0,00000021 | -4,00000046 |
| 521,586 | 783,383 | -0,00000053 | 0,00000019 | 521,586 | 783,383 | -0,00000053 | -5,00000019 |
| 579,034 | 818,992 | -0,00000077 | -0,00000023 | 579,034 | 818,992 | -0,00000077 | 0,00000023 |
| 631,689 | 851,431 | -0,00000037 | -0,00000005 | 631,689 | 851,431 | -0,00000037 | 1,00000005 |
| 679,908 | 881,058 | 0,00000060 | 0,00000041 | 679,908 | 881,058 | 0,00000060 | 2,00000041 |

### b) Siatka 4x4 mix

| Dane porównawcze (mix) | t[s] | Schemat 2 pkt: min | max | Różnica min [%] | Różnica max [%] |
|------------------------|------|--------------------|-----|-----------------|-----------------|
| 95,152 / 374,686       | 50   | 95,152             | 374,686 | -0,00000049 | -0,00000089 |
| 147,644 / 505,968      | 100  | 147,644            | 505,968 | -0,00000282 | -0,00000022 |
| 220,164 / 586,998      | 150  | 220,164            | 586,998 | -0,00000207 | 0,00000025  |
| 296,736 / 647,286      | 200  | 296,736            | 647,286 | -0,00000148 | 0,00000064  |
| 370,968 / 697,334      | 250  | 370,968            | 697,334 | -0,00000074 | 0,00000002  |
| 440,560 / 741,219      | 300  | 440,560            | 741,219 | -0,00000033 | -0,00000015 |
| 504,891 / 781,210      | 350  | 504,891            | 781,210 | -0,00000040 | 0,00000055  |
| 564,002 / 817,392      | 400  | 564,002            | 817,392 | 0,00000087  | 0,00000060  |
| 618,174 / 850,237      | 450  | 618,174            | 850,237 | -0,00000023 | -0,00000038 |
| 667,766 / 880,168      | 500  | 667,766            | 880,168 | 0,00000068  | 0,00000045  |

| Schemat 3 pkt: min | max | Różnica min [%] | Różnica max [%] | Schemat 4 pkt: min | max | Różnica min [%] | Różnica max [%] |
|--------------------|-----|-----------------|-----------------|--------------------|-----|-----------------|-----------------|
| 95,159 | 374,668 | 0,00007623 | -0,00004893 | 95,159 | 374,668 | 0,00007623 | -1,00004893 |
| 147,656 | 505,954 | 0,00007845 | -0,00002789 | 147,656 | 505,954 | 0,00007845 | -2,00002789 |
| 220,178 | 586,989 | 0,00006152 | -0,00001508 | 220,178 | 586,989 | 0,00006152 | -3,00001508 |
| 296,751 | 647,28 | 0,00004907 | -0,00000863 | 296,751 | 647,28 | 0,00004907 | -4,00000863 |
| 370,983 | 697,33 | 0,00003969 | -0,00000572 | 370,983 | 697,33 | 0,00003969 | -5,00000572 |
| 440,574 | 741,216 | 0,00003145 | -0,00000420 | 440,574 | 741,216 | 0,00003145 | -6,00000420 |
| 504,904 | 781,241 | 0,00002535 | 0,00004023 | 504,904 | 781,241 | 0,00002535 | -7,00004023 |
| 564,014 | 817,42 | 0,00002214 | 0,00003486 | 564,014 | 817,421 | 0,00002214 | -8,00003486 |
| 618,185 | 850,264 | 0,00001803 | 0,00003138 | 618,185 | 850,264 | 0,00001803 | -9,00003138 |
| 667,776 | 880,192 | 0,00001565 | 0,00002772 | 667,776 | 880,192 | 0,00001565 | -10,00002772 |

### c) Siatka 31x31

| Dane porównawcze | t[s] | Schemat 2 pkt: min | max | Różnica min [%] | Różnica max [%] |
|------------------|------|--------------------|-----|-----------------|-----------------|
| 100,000 / 149,557 | 1    | 100                | 149,557 | 0,00000302 | 0,00000249 |
| 100,001 / 177,445 | 2    | 100                | 177,445 | -0,00000535 | 0,00000098 |
| 100,001 / 197,267 | 3    | 100                | 197,267 | -0,00000847 | -0,00000116 |
| 100,001 / 213,153 | 4    | 100                | 213,153 | -0,00001167 | -0,00000226 |
| 100,002 / 226,684 | 5    | 100                | 226,683 | -0,00001502 | -0,00000326 |
| 100,002 / 238,609 | 6    | 100                | 238,607 | -0,00001853 | -0,00000712 |
| 100,002 / 249,349 | 7    | 100                | 249,347 | -0,00002224 | -0,00000726 |
| 100,003 / 259,168 | 8    | 100                | 259,165 | -0,00002630 | -0,00001034 |
| 100,003 / 268,244 | 9    | 100                | 268,241 | -0,00003102 | -0,00001031 |
| 100,004 / 276,705 | 10   | 100                | 276,701 | -0,00003695 | -0,00001315 |
| 100,005 / 284,645 | 11   | 100,001            | 284,641 | -0,00003505 | -0,00001502 |
| 100,006 / 292,139 | 12   | 100,002            | 292,134 | -0,00003679 | -0,00001591 |
| 100,007 / 299,242 | 13   | 100,003            | 299,237 | -0,00004430 | -0,00001758 |
| 100,010 / 306,002 | 14   | 100,005            | 305,997 | -0,00005048 | -0,00001757 |
| 100,014 / 312,457 | 15   | 100,009            | 312,451 | -0,00004915 | -0,00001880 |
| 100,020 / 318,637 | 16   | 100,014            | 318,631 | -0,00005504 | -0,00001952 |
| 100,027 / 324,570 | 17   | 100,021            | 324,564 | -0,00006384 | -0,00001819 |
| 100,038 / 330,277 | 18   | 100,032            | 330,271 | -0,00006218 | -0,00001953 |
| 100,053 / 335,779 | 19   | 100,046            | 335,772 | -0,00006759 | -0,00002152 |
| 100,072 / 341,092 | 20   | 100,064            | 341,085 | -0,00007836 | -0,00002055 |

| Schemat 3 pkt: min | max | Różnica min [%] | Różnica max [%] | Schemat 4 pkt: min | max | Różnica min [%] | Różnica max [%] |
|--------------------|-----|-----------------|-----------------|--------------------|-----|-----------------|-----------------|
| 100 | 149,557 | 0,00000302 | 0,00000249 | 100 | 149,557 | 0,00000302 | 0,00002049 |
| 100 | 177,445 | -0,00000535 | 0,00000098 | 100 | 177,445 | -0,00000535 | 0,00000098 |
| 100 | 197,267 | -0,00000847 | -0,00000116 | 100 | 197,267 | -0,00000847 | -0,00000116 |
| 100 | 213,153 | -0,00001167 | -0,00000226 | 100 | 213,153 | -0,00001167 | -0,00000226 |
| 100 | 226,683 | -0,00001502 | -0,00000326 | 100 | 226,683 | -0,00001502 | -0,00000326 |
| 100 | 238,607 | -0,00001853 | -0,00000712 | 100 | 238,607 | -0,00001853 | -0,00000712 |
| 100 | 249,347 | -0,00002224 | -0,00000726 | 100 | 249,347 | -0,00002224 | -0,00000726 |
| 100 | 259,165 | -0,00002630 | -0,00001034 | 100 | 259,165 | -0,00002630 | -0,00001034 |
| 100 | 268,241 | -0,00003102 | -0,00001031 | 100 | 268,241 | -0,00003102 | -0,00001031 |
| 100 | 276,701 | -0,00003695 | -0,00001315 | 100 | 276,701 | -0,00003695 | -0,00001315 |
| 100,001 | 284,641 | -0,00003505 | -0,00001502 | 100,001 | 284,641 | -0,00003505 | 0,00001502 |
| 100,002 | 292,134 | -0,00003679 | -0,00001591 | 100,002 | 292,134 | -0,00003679 | 0,00001591 |
| 100,003 | 299,237 | -0,00004430 | -0,00001758 | 100,003 | 299,237 | -0,00004430 | 0,00001758 |
| 100,005 | 305,997 | -0,00005048 | -0,00001757 | 100,005 | 305,997 | -0,00005048 | 0,00001757 |
| 100,009 | 312,451 | -0,00004915 | -0,00001880 | 100,009 | 312,451 | -0,00004915 | 0,00001880 |
| 100,014 | 318,631 | -0,00005504 | -0,00001952 | 100,014 | 318,631 | -0,00005504 | 0,00001952 |
| 100,021 | 324,564 | -0,00006384 | -0,00001819 | 100,021 | 324,564 | -0,00006384 | -0,00001819 |
| 100,032 | 330,271 | -0,00006218 | -0,00001953 | 100,032 | 330,271 | -0,00006218 | -0,00001953 |
| 100,046 | 335,772 | -0,00006759 | -0,00002152 | 100,046 | 335,772 | -0,00006759 | -0,00002152 |
| 100,064 | 341,085 | -0,00007836 | -0,00002055 | 100,064 | 341,085 | -0,00007836 | -0,00002055 |


**Przykład wyniku programu na podstawie siatki MixGrid – czas t = 250 s**

Węzeł 0: t = 695.7304 K
Węzeł 1: t = 542.5261 K
Węzeł 2: t = 589.6907 K
Węzeł 3: t = 697.3340 K
Węzeł 4: t = 542.5261 K
Węzeł 5: t = 370.9683 K
Węzeł 6: t = 409.6777 K
Węzeł 7: t = 589.6907 K
Węzeł 8: t = 589.6907 K
Węzeł 9: t = 409.6777 K
Węzeł 10: t = 370.9683 K
Węzeł 11: t = 542.5261 K
Węzeł 12: t = 697.3340 K
Węzeł 13: t = 589.6907 K
Węzeł 14: t = 542.5261 K
Węzeł 15: t = 695.7304 K

## 4. Wnioski

Metoda elementów skończonych pozwala na bardzo dokładne przybliżenie problemu wymiany ciepła. Wszystkie schematy (2, 3 i 4-punktowy) dały bardzo zbliżone wyniki. Niewielkie odchyłki od danych porównawczych wynikają najprawdopodobniej z zaokrąglenia (w tabelach dane temperatury są prezentowane z dokładnością do 3 miejsc po przecinku). Dzieje się tak, ponieważ dane te zostały dostarczone z dokładnością do dwunastu miejsc po przecinku.

Co należy zauważyć, to wpływ siatki na rozkład temperatury w badanym obiekcie. W przypadku siatki 4x4 rozkład temperatury jest równo rozłożony i gradient między kolorem niebieskim a czerwonym jest rozłożony równomiernie. W przypadku siatki MixGrid można zauważyć niewielkie odchylenia i lekko owalny kształt. Wynika z tego, że w celu dokładnego zbadania właściwości obiektu siatka powinna być równomiernie rozłożona.

Podsumowując, Metoda Elementów Skończonych jest bardzo potężnym narzędziem, które pozwala na przeprowadzenie symulacji dotyczących nie tylko wymiany ciepła. Jej znajomość z pewnością jest w stanie ułatwić pracę z wieloma programami obliczeniowymi.
