# Lipidomics (SPANISH VERSION)

# Cambios en la composición lipídica de células tiroideas bajo el control del eje TSH-CREB3L1 

### 1 Cambios en la composición lipídica bajo el estímulo con TSH
Nuestro laboratorio hipotetizó que, a partir de la estimulación de las células tiroideas a través de la TSH, al mismo tiempo que se produce una adaptación de la vía secretora debido a un incremento en la síntesis protéica, también debe darse una reprogramación lipídica que contribuye a un aumento de las membranas involucradas en el tráfico intracelular. Para llevar a cabo el estudio de esta reprogramación, evaluamos la composición lipídica de células de tiroides de rata PCCL3 estimuladas con TSH a diferentes tiempos, a través de ensayos de lipidómica utilizando una técnica basada en cromatografía líquida acoplada a un espectro de masas. En este sentido, células PCCL3 48 h luego de haber sido plaquedas en placas P100 bajo condiciones de crecimiento, fueron deprivadas de TSH por 72 h, para luego ser estimuladas con TSH a 24 h y 48 h. A modo de control del funcionamiento del esquema de deprivación y estimulación con TSH, se tomó una pequeña fracción de estas células tratadas y se analizaron por Western blot los niveles de CREB3L1, GM130 y NIS (Figura 26), en donde se observó una caída de estas tres proteínas frente a la deprivación con TSH (0 h) en comparación con células control (Ctrl, células en medio de creciemiento). Luego de la estimulación con TSH se observa un aumento marcado de CREB3L1, GM130 y NIS, observándose de forma más clara un aumento de CREB3L1 a las 24 h mientras que NIS y GM130 tienen nieveles más altos a las 48 h (Figura 1).

<img width="391" alt="image" src="https://user-images.githubusercontent.com/48334248/165776493-815a69c1-73a0-426d-862e-0de4e1a71901.png">
Figura 1. Análisis por western blot de células PCCL3 estimuladas con TSH a diferentes tiempos: Western blot de extractos celulares de células PCCL3. Células en condiciones de crecimiento (Ctrl), deprivadas de TSH por 72 h (0 h) y estimuladas con TSH por 24 h y 48 h. Etiquetas a la derecha del western blot indican la movilidad electroforética relativa del correspondiente polipéptido de NIS según su estado de glicosilación: glicosilación inmadura (~60kDa, Banda A) y la totalmente glicosilada (~100kDa, Banda B).

Ensayos de lipidómica revelaron una alta heterogeneidad en el comportamiento de los lípidos estudiados (Figura 2), en donde se puede observar un aumento de algunos grupos de lípidos como Fosfatidilcolina (PC), Fosfatidilserina (PS) y Cardiolipinas (CL) frente a la deprivación de TSH (0 h, Figura 27). Por otro lado, Diacilgliceroles (DAG) y Triacilgliceroles (TAG) son de los grupos de lípidos que más aumentan frente a la estimulación con TSH (48 h, Figura 2), particularmente luego de 48 h de estimulación. Un comportamiento similar a los DAG y TAG puede observarse en LBPA, un grupo de lípidos asociados a endosomas tardíos, y en el grupo de los esfingolípidos (SLs), que involucran a las ceramidas (CER), glucosilceramidas (GlcCer) y esfingomielinas (SM) (Figura 2).
Los esfingolípidos comprenden un pequeño grupo de lípidos que contribuyen a la función del Golgi, en parte esto se debe a su rol en señalamiento, transporte de proteínas y transformación de membranas (Halter et al., 2007). Datos obtenidos de los ensayos de lipidómica revelan en detalle los cambios observados en algunos esfingolípidos. Como se observa en la Figura 27, SM disminuye sus niveles considerablemente luego de la deprivación de TSH, para luego retomar sus niveles normales (e incluso aumentarlos) frente a una estimulación con TSH, algo muy similar ocurre en el caso de GlcCer, mientras que CER y otros SLs tienen un comportamiento más heterogéneo (Figura 2).

<img width="259" alt="image" src="https://user-images.githubusercontent.com/48334248/165777078-a7f7d434-0233-43cb-9cc4-953b3144eff9.png">

Figura 2. Ensayos de lipidómica en PCCL3 revelan un comportamiento heterogéneo en los cambios de lípidos: Células PCCL3 en condiciones de crecimiento (Ctrl), deprivadas de TSH por 72 h (0 h), y estimuladas con TSH por 24 h y 48 h. Lípidos fueron extraidos utilizando metil-tert-butil éter (MTBE) como solvente y luego introducidos en un equipo de cromatografía líquida acoplado a espectrometría de masas (LC-MS). Heatmap del dataset obtenido conteniendo todos los lípidos analizados en un experimento con 3 replicas biológicas. 

Para estudiar más en detalle la información obtenida de los ensayos de lipidómica, realizamos un análisis comparativo de la condición de deprivación de TSH por 72 h (0 h) versus la estimulación de TSH por 48 h (Figura 3). Por medio de un volcano plot identificamos diferentes grupos de lípidos que se encuentran downregulados frente a la estimulación con TSH, como es el caso de PC, LysoPC y PS (Figura 3 A). En contraposición, un grupo aún más abundante de lípidos se encuentra upregulado frente al estímulo con TSH, estos son TAG, DAG, GlcCer, LBPA, SM, CER y PE (Figura 3 A). Análisis de ontología mediante el uso del software “LION: Lipid Ontology Enrichment Analysis” (http://www.lipidontology.com) (Molenaar et al., 2019) reveló que aquellas vías metabólicas enriqucedias frente a la estimulación con TSH pertenecían a la formación de depósitos lipídicos, “lipids droplets”, triacilglicerol y diacilglicerol principalmente (Figura 3 B), lo que va en línea con reportes previos en donde se estudió el efecto de la TSH en adipocitos (Ma et al., 2015). En un segundo lugar se puede observar que la vía metabólica de hexosilceramida, precursor de varios esfingolípidos, se encuentra enriquecida frente al estímulo con TSH (Figura 3 B). Entre las vías metabólicas que se encuentran enriquecidas en los lípidos downregulados frente al estímulo con TSH predominan el de las glicerofosfocolinas y fosfatidilcolinas (Figura 3 B), los cuales son precursores de los lípidos generados frente al estímulo con TSH.

<img width="482" alt="image" src="https://user-images.githubusercontent.com/48334248/165777432-1cdb029f-80c8-4330-9ef1-fb3421708e3f.png">

Figura 3. Volcano plot y Lipid Ontology Analysis revelaron distintas vías metabólicas reguladas frente al estímulo con TSH : (A) Volcano plot utilizando el mismo dataset que en el Heatmap de Figura 27, se seleccionó un punto de corte de fold change (fc) de 4, o log2 del fc de +/- 2 para identificar aquellos lípidos cuyos cambios eran significativos. (B) Lipid Ontology Analysis de los lípidos cuyos cambios eran significativos utilizando el software LION (http://www.lipidontology.com).


