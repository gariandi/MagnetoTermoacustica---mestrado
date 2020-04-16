%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 15-04-2020  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Gabriel R. de Andrade Silva - gabriel.silva@aluno.ufabc.edu.br / unnolab777@gmail.com  %%%%%%%                                   
%%%% Mestrado em F�sica - UFABC - 2019                                                      %%%%%%%      
%%%% Investiga��o anal�tica e num�rica do efeito termoac�stico com campo magn�tico          %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Orientador: Francisco Eugenio Mendon�a da Silveira -  francisco.silveira@ufabc.edu.br  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% esta fun��o re�ne os par�metros de input para computar perfil de
%%% press�o peff(z), do efeito magneto termoac�stico

function q = geraInput(p0,L,meio,rmax,a,b,Ta,Tb,z0,desvioT,N,n_,perfil_T,B0)
%02-12-2019
%Monta structure S de input para PLOTA_ARTIGOV5_3(S)

S = [];

S.p0 = p0; 
S.L = L;
S.meio = meio ; 
S.rmax = rmax;
S.a = a;
S.b = b;
S.Ta = Ta;
S.Tb = Tb;
S.z0 = z0; 
S.desvioT = desvioT; 
S.N = N;
S.n_ = n_;
S.perfil_T = perfil_T;
S.B0 = B0;


q = S;

end