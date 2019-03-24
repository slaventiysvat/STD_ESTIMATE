function dct = DCT_ref(in_arr)

    N=length(in_arr);
    j=1;
    w=1;
    s=0;
    dct = zeros;
    for k=0:1:N-1
        x1=sqrt((2/N));
        if k==0
            x1=1/(sqrt(N));
        end
        for n=0:1:N-1
            s=s + in_arr(j)*cos(((2*n+1)*k*pi)/(2*N));
            j=j+1;
        end
        dct(w)=s*x1;
        w=w+1;
        j=1;
        s=0;
    end

end