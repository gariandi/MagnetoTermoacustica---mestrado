%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 15-04-2020  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Gabriel R. de Andrade Silva - gabriel.silva@aluno.ufabc.edu.br / unnolab777@gmail.com  %%%%%%%                                   
%%%% Mestrado em Física - UFABC - 2019                                                      %%%%%%%      
%%%% Investigação analítica e numérica do efeito termoacústico com campo magnético          %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Orientador: Francisco Eugenio Mendonça da Silveira -  francisco.silveira@ufabc.edu.br  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Este script chama a última versão dos códigos para computar o perfil da pressão
%%% manométrica efetiva ao longo numa coluna de fluído condutor, sujeito a um gradiente de temperatura e
%%% um campo magnético externo. O método é discutido em detalhe na dissertação, mas
%%% baseia-se resumidamente numa extensão do modelo de Rott (1969).

%%% Inputs do programa: a amplitude da oscilação de pressão (p0), o
%%% comprimento do tubo(L), o meio (pode ser ar ou gálio), o raio do
%%% tubo(rmax), as condições de contorno (nros binários a e b), as
%%% temperaturas nas extremidades (Ta e Tb) e alguns parâmetros referentes
%%% ao perfil de temperatura imposto externamente.

%%% Output: (ideal) gráfico peff x z . Outros gráficos foram projetados no meio
%%% do caminho, para averiguação da coerência do método. O gráfico final
%%% ainda está um pouco difícil de obter com um só script, por isso estamos
%%% abrindo a colaboração para obtenção do mesmo neste repositório público:
%%%   https://github.com/gariandi/MagnetoTermoacustica---mestrado
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Valores típicos dos parâmetros de input usados na dissertação
%(considerados em unidades do SI)
p0 = .025 ; L = 1; meio = 'Ga' ; rmax = .0254 ; a = 1; b=0; %(tubo aberto em L = 0 e fechado em L = 1;
Ta = 0; Tb=600; z0 = .25; desvioT=L/100;N=10; n_=100; perfil_T=4; B0 = 10E-5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Inicio da solução %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = geraInput(p0,L,meio,rmax,a,b,Ta,Tb,z0,desvioT,N,n_,perfil_T,B0);
PLOTA_ARTIGOV5_3_completo(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  fim da solução %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%