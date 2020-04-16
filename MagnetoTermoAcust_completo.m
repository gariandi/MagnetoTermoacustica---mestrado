%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 15-04-2020  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Gabriel R. de Andrade Silva - gabriel.silva@aluno.ufabc.edu.br / unnolab777@gmail.com  %%%%%%%                                   
%%%% Mestrado em F�sica - UFABC - 2019                                                      %%%%%%%      
%%%% Investiga��o anal�tica e num�rica do efeito termoac�stico com campo magn�tico          %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Orientador: Francisco Eugenio Mendon�a da Silveira -  francisco.silveira@ufabc.edu.br  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Este script chama a �ltima vers�o dos c�digos para computar o perfil da press�o
%%% manom�trica efetiva ao longo numa coluna de flu�do condutor, sujeito a um gradiente de temperatura e
%%% um campo magn�tico externo. O m�todo � discutido em detalhe na disserta��o, mas
%%% baseia-se resumidamente numa extens�o do modelo de Rott (1969).

%%% Inputs do programa: a amplitude da oscila��o de press�o (p0), o
%%% comprimento do tubo(L), o meio (pode ser ar ou g�lio), o raio do
%%% tubo(rmax), as condi��es de contorno (nros bin�rios a e b), as
%%% temperaturas nas extremidades (Ta e Tb) e alguns par�metros referentes
%%% ao perfil de temperatura imposto externamente.

%%% Output: (ideal) gr�fico peff x z . Outros gr�ficos foram projetados no meio
%%% do caminho, para averigua��o da coer�ncia do m�todo. O gr�fico final
%%% ainda est� um pouco dif�cil de obter com um s� script, por isso estamos
%%% abrindo a colabora��o para obten��o do mesmo neste reposit�rio p�blico:
%%%   https://github.com/gariandi/MagnetoTermoacustica---mestrado
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Valores t�picos dos par�metros de input usados na disserta��o
%(considerados em unidades do SI)
p0 = .025 ; L = 1; meio = 'Ga' ; rmax = .0254 ; a = 1; b=0; %(tubo aberto em L = 0 e fechado em L = 1;
Ta = 0; Tb=600; z0 = .25; desvioT=L/100;N=10; n_=100; perfil_T=4; B0 = 10E-5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Inicio da solu��o %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = geraInput(p0,L,meio,rmax,a,b,Ta,Tb,z0,desvioT,N,n_,perfil_T,B0);
PLOTA_ARTIGOV5_3_completo(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  fim da solu��o %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%