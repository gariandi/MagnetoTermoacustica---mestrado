%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 15-04-2020  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Gabriel R. de Andrade Silva - gabriel.silva@aluno.ufabc.edu.br / unnolab777@gmail.com  %%%%%%%                                   
%%%% Mestrado em Física - UFABC - 2019                                                      %%%%%%%      
%%%% Investigação analítica e numérica do efeito termoacústico com campo magnético          %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Orientador: Francisco Eugenio Mendonça da Silveira -  francisco.silveira@ufabc.edu.br  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Funções chamadas no script MagnetoTermoAcust_completo.m para computar o perfil da pressão
%%% manométrica efetiva ao longo numa coluna de fluído condutor, sujeito a um gradiente de temperatura e
%%% um campo magnético externo. O método é discutido em detalhe na dissertação, mas
%%% baseia-se resumidamente numa extensão do modelo de Rott (1969).


function q = PLOTA_ARTIGOV5_3_completo(S)
%02-12-2019
%V5_3: Introduzindo sugestão do prof. Zwinglio: input via struture S

%p0 = S.p0; 
L = S.L;
%meio = S.meio ; rmax = S.rmax ; a = S.a ; b = S.b;%Ta = S.Ta ; Tb = S.Tb ;
%z0 = S.z0 ; desvioT = S.desvioT; N = S.N ; n_ = S.n_; B0 = S.B0;

%06-08-2019
%Calcula correção magnética para pressão em r = 0;
%Agora com correção magnética!!!
%usando modelagem por patamares

dz = L/100;
z_vetor = 0:dz:L ;  %na realidade seria zeta ????
p = peffV4(S,z_vetor);

modulo = abs( p ) ;  fase = angle( p );

figure();
subplot(1,2,1); plot(z_vetor,modulo); grid on; xlabel('zeta'); ylabel('|p(z)|  (Pa)'); title('valor absoluto pressão manométrica');
subplot(1,2,2); plot(z_vetor,modulo); grid on; xlabel('zeta'); ylabel('arctg[Im(p)/Re(p)]  (rad)'); title('fase da pressão manométrica');

q = 1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = peffV4(S,z1)
%02-12-2019
%V4: Introduzindo sugestão do prof. Zwinglio: input via struture S
%p0 = S.p0; L = S.L; meio = S.meio ; rmax = S.rmax ; a = S.a ; b = S.b;
%Ta = S.Ta ; Tb = S.Tb ; z0 = S.z0 ; desvioT = S.desvioT; N = S.N ; n_ = S.n_;
%B0 = S.B0;

%06-08-2019
%PARA ARTIGO
% essa função determina peff(z) compondo o modelo discretizado de perfil_T

amplitudes = ResModelArtDiscretoV2(S);

size(amplitudes)

size(z1)

load Degraus; load ks;

%[i,j] = find( Degraus(:,1)<=z1 & Degraus(:,2)>z1 ) 
%isempty(i)
%size(amplitudes(i))

%p = @(z) amplitudes(i)*exp(ks(i,1)*z) + amplitudes(i+1)*exp(ks(i,2)*z) ;

%q = p(z1);

i = 1;
p = @(z) ( heaviside(z-Degraus(i,1)) - heaviside(z - Degraus(i,2)) ).*(  amplitudes(i)*exp(ks(i,1)*z) + amplitudes(i+1)*exp(ks(i,2)*z)  ) ;

for i = 2:N

    p = @(z) p(z) + ( heaviside(z-Degraus(i,1)) - heaviside(z - Degraus(i,2)) ).*(  amplitudes(i)*exp(ks(i,1)*z) + amplitudes(i+1)*exp(ks(i,2)*z)  ); 

end
q = p(z1);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = ResModelArtDiscretoV2(S)
%V2: Introduzindo sugestão do prof. Zwinglio: input via struture S
p0 = S.p0;
%L = S.L; meio = S.meio ; rmax = S.rmax ; 
a = S.a ; b = S.b;
%Ta = S.Ta ; Tb = S.Tb ; z0 = S.z0 ; desvioT = S.desvioT; 
N = S.N ; %n_ = S.n_; %B0 = S.B0;

%06-08-2019
%resolve sistema M.p = b para achar amplitudes de oscilação p 

%Discretizando o perfil de temperatura em N pedaços e resolvendo por partes
% n_ é a discretização para tirar o valor médio de cada pedaço do T0(z)
% anteriormente n era o  harmonico (agora estamos setando n = 1 para fundamental)

M = GeraMatrizMV2(S);
plota_M(M);

amplitudes = zeros(2*N,1);
b_vetor = zeros(2*N,1);

%conds de contorno possíveis (parametrizando em a e b)
if a==0 & b==1          
    b_vetor(1) = p0 ;
else
    if a==0 & b==0
        b_vetor(1) = p0 ;     b_vetor(2*N) = p0 ;
    else
        if a==1 & b==0
            b_vetor(2*N) = p0 ;
        else
            if a==1 & b==1
                b_vetor(N) = p0 ;
             else
                disp('Valores de a e/ou b incompatíveis com entrada, esperava-se zero ou 1');
            end
        end
    end
end

amplitudes = linsolve(M,b_vetor);  



q = amplitudes;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = GeraMatrizMV2(S)
%V2: Introduzindo sugestão do prof. Zwinglio: input via struture S
%p0 = S.p0; L = S.L; meio = S.meio ; rmax = S.rmax ; a = S.a ; b = S.b;
%Ta = S.Ta ; Tb = S.Tb ; z0 = S.z0 ; desvioT = S.desvioT; 
N = S.N ; %n_ = S.n_; B0 = S.B0;

%08-08-2019
% ALTERADA função ks da linha 13

%gerador de Matriz para sistema linear de amplitudes do 

%05-08-2019
%Discretizando o perfil de temperatura em N pedaços e resolvendo por partes
% n_ é a discretização para tirar o valor médio de cada pedaço do T0(z)
% anteriormente n era o  harmonico (agora estamos setando n = 1 para fundamental)


%ks = ks_T_discretizada(N,n_,a,b,L,meio,rmax,Ta,Tb,perfil_T,desvioT,z0); 
ks = ks_T_discretizadaV3(S);  %ALTERADO para V3 enm 02-12-2019 
save ks;
n_ks = size(ks); n_ks = n_ks(1);

load Degraus;
zs = Degraus(:,2) ; 
zs(length(zs)) = [] ; % assim só ficam pontos de fronteira internos (tamanho (N-1))

plota_ks(Degraus,ks);

%matriz M terá 2N linhas, amplitudes virão do sistema linear M.u = b

M = zeros(2*N,2*N) ;

M(1,1:2) = [ 1   1 ] ;

M(2*N,(2*N-1):(2*N)) = [ exp(ks(n_ks,1)*L)    exp(ks(n_ks,2)*L) ]  ;

for bloco = 1:(N-1)
    
    z = zs(bloco);
    
    ka_mais = ks(bloco,1) ; ka_menos = ks(bloco,2);
    kb_mais = ks(bloco+1,1) ; kb_menos = ks(bloco+1,2) ; 

    imin = 2*bloco  ;  imax = 2*bloco+1 ;
    jmin = (2*bloco-1) ; jmax = (2*bloco+2) ;
       
    M(imin:imax,jmin:jmax) =  [     exp(ka_mais*z)          exp(ka_menos*z)              -exp(kb_mais*z)           -exp(kb_menos*z)   ;
                                ka_mais*exp(ka_mais*z)   ka_menos*exp(ka_menos*z)   -kb_mais*exp(kb_mais*z)   -kb_menos*exp(kb_menos*z)   ];
                            
                            
end

q = M;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = plota_M(M)

figure();

%subplot(1,2,1); surf(real(M)); view(0,270); colorbar;
subplot(1,2,1); surf(real(M)); view(0,-90); colorbar; %alterado para rodar no octave
title('parte real de M') ; 

%subplot(1,2,2); surf(imag(M)); view(0,270); colorbar;
subplot(1,2,2); surf(imag(M)); view(0,-90); colorbar; %alterado para rodar no octave
title('parte imaginária de M') ;

q = 1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = ks_T_discretizadaV3(S)
%V3: Introduzindo sugestão do prof. Zwinglio: input via struture S
p0 = S.p0; L = S.L; meio = S.meio ; rmax = S.rmax ; a = S.a ; b = S.b;
Ta = S.Ta ; Tb = S.Tb ; z0 = S.z0 ; desvioT = S.desvioT; N = S.N ; n_ = S.n_;
perfil_T = S.perfil_T; B0 = S.B0;

%08-08-2019
% setando ABC da linha 32 para usar coefs_artigoV3 (perfil exclusivamente cosh) 
n = 1; % setando harmônico fundamental


%08-08-2019
% V2: calcula k_mais e k_menos por Bhaskara, sem usar função solve()

%05-08-2018: linha 11 atualizada com Discretiza_T_V2
% n = harmonico (1 para fundamental)
% n_ = nro divisoes para fazer média dos intervalos

%17-07-2019: procedendo solução por discretização do perfil de temperatura,
% e divisão do dominio em sessões

% considera perfil dado por T0_(z) salva em outro arquivo

ambiente = parametros_meio(meio);
gama_cte = ambiente(3);
%save Ta; save Tb; save perfil_T; save desvioT; save z0; save gama_cte;

Degraus = Discretiza_T_V3(L,N,n_,Ta,Tb,perfil_T,desvioT,gama_cte,z0);   %V2 e V3 são equivalentes
save Degraus;

dz = L/n_;

ks = zeros(N,2);

 for m = 1:N
    
    z = Degraus(m,1):dz:Degraus(m,2) ;
  %  ABC = coefs_artigo(a,b,L,meio,rmax,1,Ta,Tb,perfil_T,desvioT,z0,z) ; %setando harmonico fundamental ( n = 1 ), para não dar confusão com outro n
  %  ABC = coefs_artigoV3(a,b,L,meio,rmax,n,Ta,Tb,z);
     ABC = coefs_artigoV6(a,b,L,meio,rmax,Ta,Tb,z);
  
    A = mean(ABC(:,1)) ;  B = mean(ABC(:,2)) ;  C = mean(ABC(:,3)) ;
        
    k_mais = ( -B + sqrt(B.^2 - 4*A*C) )/(2*A) ;
    k_menos = ( -B - sqrt(B.^2 - 4*A*C) )/(2*A) ;   
    k = [k_mais k_menos];
    
    ks(m,:) = k ;
    
 end


q = ks;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = coefs_artigoV6(a,b,L,meio,rmax,Ta,Tb,z1)
%21-08-2019
%V6: trocando por besseli e retirando unidade imaginária do
%argumento

%11-08-2019
%V5: rastreando todas funções via plot - a procura do NaN !!!

%08-08-2019
%V4: redefinindo coeficientes SEM SINGULARIDADES

n=1;

%V3: define T(z) por perfil_temperaturaV2() - exclusivamente cosh !!!

%21-06-2019: esta função devolve os três coeficientes A, B e C da equação
%em um vetor linha, calculados para z = z1 .

%%%%%%%%%parâmetros do meio%%%%%%%%%%%%%%%%%
Ambiente = parametros_meio(meio); % [  R   patm   gama_cte  M   K_ni  beta_ni  Pr  ] ;
R = Ambiente(1) ;  patm = Ambiente(2) ;  gama_cte = Ambiente(3) ;  
M = Ambiente(4); 
K_ni = Ambiente(5) ;  beta_ni = Ambiente(6) ;  Pr = Ambiente(7) ;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save R; save M; save patm; save K_ni; save beta_ni; save Pr; save gama_cte;
%load dz;
dz = L/100000;

% definindo perfil T0(z)
T = @(z) ( perfil_temperaturaV2(Ta,Tb,gama_cte,z) );

gradT = @(z) ( (T(z+dz/2)-T(z-dz/2))/dz );
theta = @(z) ( (T(z).^-1).*gradT(z) );

lambda = lamb(a,b,L,n);
k = 2*pi/lambda; 

vsom = @(z) ( sqrt(gama_cte*R.*T(z)./M) );
omega = @(z) ( k*vsom(z) );

ni = @(T) ( K_ni*T.^beta_ni );  %CHECAR VALORES DE K_ni e beta_ni para cada meio

%%%%acrescentado aqui para simular artigo em %%%%11-04-2019%%%%%%%%%%%%%%%%

figure();
z_vetor = 0:dz:L;
lc = @(z) ( sqrt( omega(z)./ni(T(z)) ) );  %INVERSO da camada limite
subplot(2,3,1); plot(z_vetor,lc(z_vetor)); grid on; title('inverso da camada limite');

qsi0 = @(z) ( sqrt(i*Pr)*lc(z)*rmax );   %qsi0 = @(z) rmax*sqrt(i*Pr*omega(z)/ni(z));
subplot(2,3,2); plot(z_vetor,real(1./qsi0(z_vetor))); grid on; title('1/qsi0(z)'); hold on;
                          plot(z_vetor,imag(1./qsi0(z_vetor))); legend('parte real','parte imag');
qsi1 = @(z) ( sqrt(i)*lc(z)*rmax );
subplot(2,3,3); plot(z_vetor,1./qsi1(z_vetor)); grid on; title('1/qsi1(z)');

subplot(2,3,4); plot(z_vetor,real(besseli(1,qsi0(z_vetor)))); grid on; hold on; title('funções de bessel J1 e 1/J0');
                          plot(z_vetor,imag(besseli(1,qsi0(z_vetor)))); hold on;
                          plot(z_vetor,real(1./besseli(0,qsi0(z_vetor)))); hold on;
                          plot(z_vetor,imag(1./besseli(0,qsi0(z_vetor)))); legend('Re[J1(i.qsi0)]','Im[J1(i.qsi0)]','Re[1/J0(i*qsi0)]','Im[1/J0(i*qsi0)]');
                          
fi_sigma = @(z) ( 1 + (2./qsi0(z)).*besseli(1,qsi0(z))./besseli(0,qsi0(z)) );
subplot(2,3,5); plot(z_vetor,real(fi_sigma(z_vetor))); grid on; hold on; title('fi_sigma(z)');
                          plot(z_vetor,imag(fi_sigma(z_vetor)));  legend('parte real','parte imag');


fi_1 = @(z) ( 1 + (2./qsi1(z)).*besseli(1,qsi1(z))./besseli(0,qsi1(z)) );
gradfi_1 = @(z) ( (fi_1(z+dz/2)-fi_1(z-dz/2))/dz );
subplot(2,3,6); plot(z_vetor,real(fi_sigma(z_vetor))); grid on; hold on; title('fi_1(z)');
                          plot(z_vetor,imag(fi_sigma(z_vetor)));  legend('parte real','parte imag'); hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%VAMOS MULTIPLICAR TODOS COEFICIENTES PELAS FUNÇÕES NO
%DENOMINADOR:¨(  qsi0(z).*besselj(0,i*qsi0(z)).*besselj(0,i*qsi1(z))  ).*

%%%REdefinindo coeficientes da equação
% d/dz[A(z).dp/dz] - B(z).dp/dz + C(z).p = 0
%vsom = @(z) sqrt(gama_cte*R.*T(z)./M);
%A = @(z) ( (  qsi0(z).*besselj(0,i*qsi0(z)).*besselj(0,i*qsi1(z))  ).*( fi_1(z) ) );  %ao invés de %A = @(z) ((vsom(z)./omega)^2)*( 1 - f(eta0(T(z))) );
A = @(z) (   fi_1(z)  );  

%B = @(z) ( (  qsi0(z).*besselj(0,i*qsi0(z)).*besselj(0,i*qsi1(z))  ).*( gradfi_1(z) + inv(1-Pr)*( fi_sigma(z) - Pr*fi_1(z) ).*theta(z) ) ); %ao invés de %B = @(z) (vsom(z)/omega)^2 * inv(1-Pr) * ( f(sqrt(Pr)*eta0(T(z))) - f(eta0(T(z))) );
B = @(z) (  gradfi_1(z) + inv(1-Pr)*( fi_sigma(z) - Pr*fi_1(z) ).*theta(z)  ); 

%C = @(z) ( (  qsi0(z).*besselj(0,i*qsi0(z)).*besselj(0,i*qsi1(z))  ).*( gama_cte - (gama_cte - 1)*fi_sigma(z) ) ); %ao invés de %C = @(z) 1 + (gama_cte-1)*f(sqrt(sigma)*eta0(T(z)));
C = @(z) (  gama_cte - (gama_cte - 1)*fi_sigma(z)  ); 
%%%%

Saida = [ A(z1)' B(z1)' C(z1)'];

q = conj(Saida);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = perfil_temperaturaV2(Ta,Tb,gama_cte,z)
%08-08-2019
%V2: perfil EXCLUSIVAMENTE cosh


%convertendo para Kelvin
Ta = Ta + 273.15 ; Tb = Tb + 273.15;

%determinando maior e menor T
if Ta>Tb
    Tmax = Ta ;  T0 = Tb;
else
    if Tb>Ta
        T0 = Ta ;   Tmax = Tb;
    else
        disp('Não existe gradiente de T. O fenômeno não pode ocorrer');
    end
end

%Perfil de Temperatura ao longo de z
  T = @(z) ( Tmax*cosh(sqrt(gama_cte-1)*z).^((2*gama_cte-1)/(gama_cte-1)) );


q = T(z);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = lamb(a,b,L,n)

disp('O comprimento do tubo correspondente é '); disp(L); disp('          metros');

if a==0 && b==0
    lambda = (2/(2*n+1))*L;
else
    if a==0 && b==1
        lambda = (4/(2*n+1))*2*L;
    else
        if a==1 && b==0
           lambda = (4/(2*n+1))*2*L;
        else
            if a==1 && b==1
                lambda = (2/n)*L;
            else
                disp('Valor(es) incompativel(is) de a e/ou b')
            end
        end
    end
end

disp('O comprimento da onda nesse tubo é aproximadamente'); disp(lambda); disp('          metros');

q = lambda;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = plota_ks(Degraus,ks)

zeta_med = mean(Degraus(:,1:2)');

figure();

subplot(1,2,1); plot(zeta_med,real(ks(:,1)),'r'); hold on ; plot(zeta_med,real(ks(:,2)),'b');
title('parte real de k'); grid on; legend('k+','k-') ;

subplot(1,2,2); plot(zeta_med,imag(ks(:,1)),'r'); hold on; plot(zeta_med,imag(ks(:,2)),'b');
title('parte imaginária de k'); grid on; legend('k+','k-') ;

q = 1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = parametros_meio(meio)

saida = zeros(1,7);

patm = 101325; %[Pa] ;
gama_cte = 1.4 ; %adimensional


switch meio
    case 'ar'

         M = 0.028964; %[kg/mol] - ar seco;          
         R = 8.3144621; %[SI]; 
        %perfil da visc cinemática X T 
         K_ni = 7*10^-10;  beta_ni = 1.7558; %(parametros do ar)
        %número de Prandtl
         Pr = 0.7 ;
         
    case 'Ga'

        M = 69.723*10^-3; %[kg/mol] - gálio
        R = (1.8E6)*M; %[SI] - ajuste arquivo "Análise Dados equação de estado do gálio.xlsx" 
        % y = 0,0003x-0,842 %ajuste excel para o gálio
        K_ni = 3E-4;  beta_ni = -0.842; %para o GÁLIO LÍQUIDO
        Pr = 0.025 ;  % GÁLIO LÍQUIDO
    
    otherwise
      
       disp('Meio não cadastrado');
        M = [] ;
        K_ni = [];  
        beta_ni = []; 
        Pr = [] ; 
        
end


%M

% saida = [  R   patm   gama_cte  M   K_ni  beta_ni  Pr  ] ;

saida = horzcat(  R ,  patm  , gama_cte , M ,  K_ni , beta_ni , Pr  ) ;

q = saida;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = Discretiza_T_V3(L,N,n,Ta,Tb,perfil_T,desvioT,gama_cte,z0)
%01-08-2019
%esta é uma versão simplificada, discretizando a passos iguais no eixo x

%31-07-2019
%agora não precisa f de entrada, setei T(z) aqui

%load Ta; load Tb; load perfil_T; load desvioT; load gama_cte; load z0;

f = @(z) perfil_temperatura(Ta,Tb,perfil_T,desvioT,gama_cte,z0,z);


% Discretiza_f: 11-07-2019
% f: função a discretizar
% L: extremo máximo do domínio ( é [0,L] )
% N: número de níveis
% n: passo de integração

dx = L/n;
x = 0:dx:L;

DX = L/N;
X = 0:DX:L; 

fmed = zeros(1,length(X)-1);

for k = 1:N
    fmed(k) = mean( f(X(k)):dx:f(X(k+1)) );
end

X_ = []; X1 = X_;

for k = 1:N
    X1 = [X(k) X(k+1)];
    X_ = vertcat(X_,X1);
end

plot(x,f(x)); grid on; hold on;
plot(X_(:,1),fmed','*'); hold on;
plot(X_(:,2),fmed','*'); hold on;
title('perfil de temperatura e discretização utilizada')
xlabel('zeta'); ylabel('T0(zeta)');

q = [X_ fmed'];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = perfil_temperatura(Ta,Tb,perfil_T,desvioT,gama_cte,z0,z)

%convertendo para Kelvin
Ta = Ta + 273.15 ; Tb = Tb + 273.15;

%determinando maior e menor T
if Ta>Tb
    Tmax = Ta ;  T0 = Tb;
else
    if Tb>Ta
        T0 = Ta ;   Tmax = Tb;
    else
        disp('Não existe gradiente de T. O fenômeno não pode ocorrer');
    end
end

%Perfil de Temperatura ao longo de z
if perfil_T==1
    T = @(z) T0 + 0.5*Tmax*( 1 + erf( (z-z0)/(sqrt(2)*desvioT) ) );  %perfil gaussiano em gradT
else if perfil_T==2
         T = @(z)  T0 + Tmax./( 1 + exp(-(z - z0)/desvioT ) ) ;  %perfil logístico de T(z)
    else if perfil_T==3
            T = @(z)  (T0 + 0.5*Tmax) + (1/pi)*Tmax*atan(2*pi*(z-z0)/desvioT);  %lorentziana em gradT
        else if perfil_T==4
               T = @(z) Tmax*cosh(sqrt(gama_cte-1)*z).^((2*gama_cte-1)/(gama_cte-1));
            else   
               error('Perfil de T(z) não cadastrado!');
               T = @(z) NaN ;
            end
        end
    end
end

save T0; save Tmax; save desvioT; save perfil_T;

q = T(z);

end
