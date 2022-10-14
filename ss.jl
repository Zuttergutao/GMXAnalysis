using Luxor

lines=readlines(ARGS[1])
ss_dict =Dict('H'=>"H", 'G'=>"H", 'I'=>"H",
           'B'=>"E", 'E'=>"E", 'P'=>"H",
           'T'=>"T", 'S'=>"T")
for (m,n) in enumerate(lines)
    if occursin("  #  RESIDUE AA",n)
        idx=m+1
        global data=lines[idx:end]
        break
    end
end
SS=[]
res=[]
for (idx,dat) in enumerate(data)
    if dat[17]==' '
        append!(SS,"T")
    else
        append!(SS,ss_dict[dat[17]])
    end
    append!(res,dat[14])
end

ss_block=[]
prev_ss=' '
for (idx,com) in enumerate(SS)
   global prev_ss
   reduced_elem=com[1]
   if idx%50!=1
      if reduced_elem!=prev_ss
         append!(ss_block,reduced_elem)
         append!(ss_block,idx)
         append!(ss_block,idx)
         prev_ss=reduced_elem
      end
      ss_block[end,end]=idx
  else
      append!(ss_block,reduced_elem)
      append!(ss_block,idx)
      append!(ss_block,(idx%50+1)*50)
      prev_ss = reduced_elem
      ss_block[end,end]=idx
  end
end
ss_block=reshape(ss_block,3,:)
tmp=Array{Any}(undef,size(ss_block,2),3)
for i in axes(ss_block,2)
   tmp[i,1]=ss_block[1,i] 
   tmp[i,2]=ss_block[2,i] 
   tmp[i,3]=ss_block[3,i]
end
ss_block=tmp

fig_w=800
fig_h=80
reslen=length(SS)
fignum=ceil(reslen/50)
t_fw=fig_w
t_fh=fig_h*fignum
coil_height=1*fig_h/50
sheet_height=7.5*fig_h/50
fs=16

Drawing(t_fw+20,t_fh,"ss.png")
background("white")
fontface("JuliaMono-ExtraBold")
for m=1:size(ss_block,1)
    ss_type,start,stop=ss_block[m,:]
    # coil
    draw_start=(fig_w*([(start%50)==0 ? 50 : (start%50)])/50)[1]
    global draw_end=(([(start%50)==0 ? (start/50)-1 : floor(start/50)])[1])*fig_h
    if ss_type=='T'
        setcolor("gray")
        Luxor.rect(Point(draw_start,(fig_h*(37.5/50)-coil_height)+draw_end), ((stop-start+1)/50)*fig_w,coil_height*2, action = :fill)
        sethue("black")
        fontsize(fs)
        text(join(res[start:stop]," "), Point(draw_start,fig_h*(25/50)+draw_end))
        fontsize(fs*0.65)
        text("|", Point(draw_start,fig_h*(15/50)+draw_end))
        text("$start", Point(draw_start,fig_h*(8/50)+draw_end))
    # helix
    elseif ss_type=='H'
        setcolor("orange")
        n_turns=ceil((stop-start+1)/1.5)
        start_t=draw_start
        wave = [Point(start_t+(((stop-start+1)*fig_w/50)/(n_turns*π))*y,fig_h*(37.5/50)+draw_end+sheet_height*sin(y)) for y in 0:π/20:(n_turns)*π]
        poly(wave, action = :stroke)
        sethue("black")
        fontsize(fs)
        text(join(res[start:stop]," "), Point(draw_start,fig_h*(25/50)+draw_end))
        fontsize(fs*0.65)
        text("|", Point(draw_start,fig_h*(15/50)+draw_end))
        text("$start", Point(draw_start,fig_h*(8/50)+draw_end))
        # sheet
    else ss_type=='E'
        setcolor("lightseagreen")
        Luxor.rect(Point(draw_start,(fig_h*(37.5/50)-sheet_height)+draw_end), ((stop-start+1)/50)*fig_w,sheet_height*2, action = :fill)
        sethue("black")
        fontsize(fs)
        text(join(res[start:stop]," "), Point(draw_start,fig_h*(25/50)+draw_end))
        fontsize(fs*0.65)
        text("|", Point(draw_start,fig_h*(15/50)+draw_end))
        text("$start", Point(draw_start,fig_h*(8/50)+draw_end))
    end
    
end
text("|", Point((fig_w*([(ss_block[end,end]%50)==0 ? 50 : ((ss_block[end,end])%50)])/50)[1],fig_h*(15/50)+draw_end))
text("$(ss_block[end,end])", Point((fig_w*([(ss_block[end,end]%50)==0 ? 50 : ((ss_block[end,end])%50)])/50)[1],fig_h*(8/50)+draw_end))
fontsize(fs)
text("C-loop", Point(5+(fig_w*([(ss_block[end,end]%50)==0 ? 50 : ((ss_block[end,end]+1)%50)])/50)[1],fig_h*(41/50)+draw_end))
finish()
