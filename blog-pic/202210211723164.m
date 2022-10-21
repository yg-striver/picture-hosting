clc;
clear ;
close ;

%% �����ļ�
load('eeg_9.mat');

%% �������
fs=128;%������
T=10;%ȡ����ʱ�䳤��
epoch=T*fs;
overlap=0.8*epoch;%����������ȡ75%����Ҫ����0.8


%% �������
r=0.15;%��ֵ,fast_ApEn��sampen

%% α����
d=10;
D=d*fs;% number of FFT points.frequency resolution=0.5
warning off;
%% ��ȡ����
for No=1%��1~58,60�����,59~100,70�������%length(TianjingData)
    EEGch1=a.pEEG1(1,:);
    %     bis=TianjingData(No).BIS;
    %     name=TianjingData(No).name;
    %     time=double(TianjingData(No).Time);
    
    length_Raweeg1=length(EEGch1);
    N=floor((length_Raweeg1-epoch)/(epoch-overlap))+1;%����ָ�곤��
    
    %�����ڴ�
    %     BIS=resample(bis,N,length(bis));%��bis������
    %     Time=resample(time,N,length(time));%��bis������
    pe=zeros(N,1);
    apen=zeros(N,1);
    sampen=zeros(N,1);
    SFS=zeros(N,1);
    Delta=zeros(N,1);
    Theta=zeros(N,1);
    Alpha=zeros(N,1);
    Beta=zeros(N,1);
    Gamma=zeros(N,1);
    RE=zeros(N,1);
    SE=zeros(N,1);
    MPF=zeros(N,1);
    SEF95=zeros(N,1);
    plzc=zeros(N,1);
    dfa=zeros(N,1);
    %     Physio=zeros(N,4);
    %% ����ָ��
    for p=1:N
        
        dataSegments=EEGch1((p-1)*(epoch-overlap)+1:(p-1)*(epoch-overlap)+epoch);
        
        pe(p,1)=perm(dataSegments,6,1);%���ݾ��飬fs=10,epoch=10*fs��m=6,tao=1
        %         apen(p,1)=fast_ApEn(dataSegments,r);%��2EEG entropy����(P4)�� from liang:N=1000,r=0.1~0.25,m=2~3,fast_ApEn equal to ApEn(when m=2)
        %         sampen(p,1)=fast_SampEn(dataSegments,2,r);%��˧:dataSegments_N=10*100,m=2,r=0.15;��2EEG entropy����(P5)��:����ApEn�ĸĽ���ȡֵ��ͬ��ApEn
        [SFS(p,1),Delta(p,1),Theta(p,1),Alpha(p,1),Beta(p,1),Gamma(p,1)] = sfsAndPower(dataSegments,epoch,T);%������Tȷ��
        %         [SEF95(p,1),MPF(p,1)] = sef95AndMPF(dataSegments,epoch,T);%������Tȷ��
        %         [RE(p,1),SE(p,1)] = SEandRE(dataSegments,T);%������Tȷ��
        plzc(p,1)=PLZC(dataSegments,4,1);%��˧:dataSegments_N=1000��10*100��,m=4,tao=1
        %         dfa(p,1)=DFA_a(dataSegments);%DFA_a����DFA�������a��a=0.5������أ�aΪ(0,0.5)�����,(0.5,1)����أ�a=1�ź�Ϊ����
        %         Physio(p,:)=Phy(No,:);%�������
        %         BIS(p,1)=bis(p);%�ز������BIS
        %         Time(p,1)=Time(p);%�ز������BIS
        p
    end
    sampleSet{No,1}=[pe,SFS,Delta,Theta,Alpha,Beta,Gamma,plzc];
    save sampleSet_pEEG1 sampleSet%ÿ����1���ͱ���һ��
%     clear pe apen sampen SFS Delta Theta Alpha Beta Gamma SEF95 MPF RE SE plzc dfa
end
