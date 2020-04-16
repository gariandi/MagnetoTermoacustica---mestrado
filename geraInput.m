%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 15-04-2020  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Gabriel R. de Andrade Silva - gabriel.silva@aluno.ufabc.edu.br / unnolab777@gmail.com  %%%%%%%                                   
%%%% Mestrado em Física - UFABC - 2019                                                      %%%%%%%      
%%%% Investigação analítica e numérica do efeito termoacústico com campo magnético          %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Orientador: Francisco Eugenio Mendonça da Silveira -  francisco.silveira@ufabc.edu.br  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% esta função reúne os parâmetros de input para computar perfil de
%%% pressão peff(z), do efeito magneto termoacústico

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