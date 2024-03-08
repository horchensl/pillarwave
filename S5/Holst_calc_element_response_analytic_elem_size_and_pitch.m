clc
clearvars
close all

Hz2MHz = 1e-6;
m2um = 1e6;

do_plot=1;
f=7e6; 

DIR_OUT = pwd;
a_var=40*1e-6;
gap_width_var=[0 20]*1e-6;

for kk=length(a_var):-1:1
    
    for ll=1:length(gap_width_var)
        % element size
        a=a_var(kk);
        b=a;
        
        px=(a+gap_width_var(ll));
        py=px;
        
        Nx=128;
        Ny=100;
      
        for ii=length(f):-1:1
            
            w=2*pi*f(ii);
            c=1500;
            
            % wavenumber axis
            kmax=1e6; Nk=1e3;
            kx=linspace(-kmax,kmax,Nk);
            ky=kx;
            
            [KX,KY]=meshgrid(kx,kx);
            
            W=sqrt(KX.^2+KY.^2)<=(w./c); % pass-band of propagation operator
            
            % calculate single element response
            A=4./(a.*b).*sin(KX.*a./2)./KX.*sin(KY.*b./2)./KY;
            
            ang=linspace(0,2*pi,1e3);
           
            % calculate full array response
            X=exp(1j.*KX.*px);
            Y=exp(1j.*KY.*py);
            
            % shift operator response
            H=(1-X.^Nx)./(1-X).*(1-Y.^Ny)./(1-Y); % array response
            
            % calculate total response
            P=A.*H;
            if do_plot
                figure(ll); clf
                imagesc(kx./m2um,ky./m2um,20*log10(abs(P)./max(abs(P(:)))),[-60 0]);
                hold on;
                h1=plot(w./c.*cos(ang)./m2um,w./c.*sin(ang)./m2um,'k');
                set(h1,'linewidth',2);
                xlabel(['kx [\mum^{-1}], f = ',num2str(f(ii)*Hz2MHz),' MHz']);
                ylabel('ky [\mum^{-1}]');
                title([{'array element response'},{['single element response, element size: ',num2str(a*m2um),' x ',num2str(b*m2um),' \mum']},{['array pitch: ',num2str(px*m2um),' x ',num2str(py*m2um),' \mum, N_x x N_y: ',num2str(Nx),' x ',num2str(Ny)]}]);
                colormap(jet)
                h1=colorbar; xlabel(h1,'[dB]');
                drawnow;
                axis image
%                 saveas(gcf,[DIR_OUT, filesep 'array_response_f_',num2str(f(ii)*Hz2MHz),'_MHz_a_',num2str(a*m2um),'_um_b_',num2str(b*m2um),'_um_pitch_',num2str([px,py]*m2um),'_Nele_',num2str([Nx Ny]),'.png'],'png');
            end
            
        end
        
    end
end



