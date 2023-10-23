function [output] = isRealizable(m_)

Hankel_Hadamard_0=zeros(length(m_)/2);
for i_=1:length(m_)/2
    for j_=1:length(m_)/2
        Hankel_Hadamard_0(i_,j_)=m_(i_+j_-1);
    end
end

Hankel_Hadamard_1=zeros(length(m_)/2);
for j_=1:length(m_)/2
    for i_=2:length(m_)/2
        Hankel_Hadamard_1(i_,j_)=m_(i_+j_);
    end
end

if det(Hankel_Hadamard_0)>=0 && det(Hankel_Hadamard_1)>=0
    output=true;
elseif det(Hankel_Hadamard_0)<0 || det(Hankel_Hadamard_1)<=0
    output=false;
end

end