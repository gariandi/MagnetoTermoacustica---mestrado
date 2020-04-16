# MagnetoTermoacustica---mestrado
Códigos matlab para simulação do efeito termoacústico e do magneto-termoacústico - modelo de Rott extendido, desenvolvido por Gabriel R. A. Silva e Francisco Eugênio M. S. , na UFABC (2016-2019), com suporte da Capes.

O script MagnetoTermoAcust_completo.m chama nossa última versão dos códigos para computar o perfil de pressão manométrica efetiva numa coluna de fluído condutor sujeito a um gradiente de temperatura e um campo magnético externo.
O problema, sua modelagem e o método de solução implementados estão discutidos resumidamente na [apresentação](https://github.com/gariandi/MagnetoTermoacustica---mestrado/blob/master/Apresentacao.pdf) e em detalhe na [dissertação](https://github.com/gariandi/MagnetoTermoacustica---mestrado/blob/master/disserta%C3%A7%C3%A3o%20Gabriel%20R%20A%20Silva%20-%20Magneto%20Termoac%C3%BAstica.pdf),
mas baseia-se resumidamente numa extensão do modelo de Rott (1969) para fenômenos termoacústicos.

Os arquivos geraInput.m e PLOTA_ARTIGOV5_3_completo.m precisam estar juntos no mesmo diretório do MagnetoTermoAcust_completo.m . Para o script funcionar basta chamar por MagnetoTermoAcust_completo no prompt do Matlab.

Inputs do programa: a amplitude da oscilação de pressão (p0), o comprimento do tubo(L), o meio (pode ser ar ou gálio), o raio do tubo(rmax), as condições de contorno (nros binários a e b), as
temperaturas nas extremidades (Ta e Tb) e alguns parâmetros referentes ao perfil de temperatura imposto externamente.

Output: (ideal) gráfico peff x z . Outros gráficos foram projetados no meio caminho, para averiguação da coerência do método. O gráfico final
ainda está um pouco difícil de obter com um só script.

A seguir uma visualização da estrutura deste projeto, com a relação entre as funções escritas para resolver o problema:

![ESTRUTURA DO PROJETO](https://github.com/gariandi/MagnetoTermoacustica---mestrado/blob/master/ESTRUTURA%20DO%20PROJETO.png)

Conforme relatado na dissertação, conseguimos obter o perfil a seguir, iterando o método de solução manualmente no prompt de comando:

![Output manual](https://github.com/gariandi/MagnetoTermoacustica---mestrado/blob/master/Output%20manual.png)

Notamos indícios de problemas com essa solução, por causa das descontinuidades. Ainda estamos por decidir se é uma questão mais computacional ou da matemática usada na construção do algoritmo.

Anexamos também um [detalhamento dos últimos bugs detectados na parte final do trabalho](https://github.com/gariandi/MagnetoTermoacustica---mestrado/blob/master/%C3%BAltimos%20bugs%20consertados.pdf) .

Plots de resultados intermediários, como o dessa matriz abaixo (sistema linear continuidade/conds contorno), parecem estar saindo normalmente, em sua maioria, nesta versão.

![resultado intermediário - matriz M](https://github.com/gariandi/MagnetoTermoacustica---mestrado/blob/master/resultado%20intermedi%C3%A1rio%20-%20matriz%20M.png)
