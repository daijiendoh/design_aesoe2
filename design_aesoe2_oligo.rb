require 'csv'
require 'bio'
require 'oj'


def calc_g60(oseq,gibbs)
	arr_gv=Array.new
	# p "oseq #{oseq}"
	slen= oseq.length
	gv=0.0
	# p oseq
	(0..slen-2).each{|i|
		gv=gv+gibbs[60.to_s][oseq[i,2].upcase]
	}
	return gv
end

def gibbs_end(oligo,gibbs,jointgibbs)
	deltagibbs=Hash.new
	(14..19).to_a.each{|elen|

		eoligo=oligo[-elen..-1]
		# p "oligo #{oligo} eoligo #{eoligo}"
		deltagibbs[elen]=calc_g60(eoligo,gibbs)-jointgibbs
		# p "#{elen} #{calc_g60(eoligo,gibbs)} #{jointgibbs} #{calc_g60(eoligo,gibbs)-jointgibbs}"
	}
	jlength= deltagibbs.sort_by{|k,v| v*v}[0][0].to_i
	return jlength
end

def gibbs_top(oligo,gibbs,jointgibbs)
	deltagibbs=Hash.new
	(14..19).to_a.each{|elen|
		eoligo=oligo[0..elen-1]
		deltagibbs[elen]=calc_g60(eoligo,gibbs)-jointgibbs
		# p "#{elen} #{calc_g60(eoligo,gibbs)} #{jointgibbs} #{calc_g60(eoligo,gibbs)-jointgibbs}"
	}
	jlength= deltagibbs.sort_by{|k,v| v*v}[0][0].to_i
	return jlength
end

jointgibbs=-13.6.to_f
# syntempl=Oj.load_file("assyntempl.oj",:mode=>:compat)

tgtfasta=ARGV[0]
syntempl=Hash.new
ff=Bio::FlatFile.auto(tgtfasta)
ff.each_entry do |f|
	acc=f.definition.chomp
	syntempl[f.definition]=f.seq.upcase
end

gibbs=Oj.load_file("gibbs.oj",:mode=>:compat)

e1range=Hash.new
syntempl.each{|acc,seq|
	slen=seq.length
	e1range[acc]=Hash.new
	i=1
	st=1
	
	olen=70
	ed=olen
	rlen=slen-st+1
	# p "#{acc} #{seq.length} #{templ_grp[acc.to_s]}" 
	# p "#{st} #{st+(ed-1)}"
	# p "slen #{slen} i #{i} olen #{olen} ed #{ed} st #{st} rlen #{rlen}"
	while rlen > 70 do
		oligo=seq[st-1..st+olen-2]
		# p oligo.length
		jendlng=gibbs_end(oligo,gibbs,jointgibbs)

		# p "#{st} #{st+(ed-1)} #{st+(ed-1)-st} #{jendlng}"
		e1range[acc][i]=[st,st+(olen-1)]
		st=st+(olen-1)-jendlng
		ed=st+olen-1
		rlen=slen-st+1
		p "#{i} #{rlen}"
		# p "slen #{slen} i #{i} olen #{olen} jendlng #{jendlng} st #{st} rlen #{rlen}"
		i+=1
		
	end

	olen=rlen	
	# p oligo.length
	jendlng=gibbs_end(oligo,gibbs,jointgibbs)
	# st=st+(olen-1)-jendlng
	oligo=seq[st-1..slen]
	# p "#{st} #{st+(ed-1)} #{st+(ed-1)-st} #{jendlng}"
	e1range[acc][i]=[st-1,slen]
	oend=st+(olen-1)
	ed=st+olen-1
	
	# p "slen #{slen} i #{i} olen #{olen} jendlng #{jendlng} st #{st} oend #{oend}"
}
# p e1range


aesoetempl=Hash.new
syntempl.each{|acc,seq|
	aesoetempl[acc]=Hash.new
	rgser=e1range[acc]
	rgser.each{|i,rg|
		aesoetempl[acc][i]=seq[rg[0]-1..rg[1]-1]

	}
	
}

File.open("aeso_f_olgo_for_check.fasta","w") do |file|
	aesoetempl.each{|acc,d1|
		d1.each{|i,seq|
			file.puts(">#{acc}_#{i}_#{seq.length}")
			file.puts(seq)
		}
	}
end

d7aeoligo=Hash.new
easet=1
aesoetempl.each{|acc,d1|
	d7aeoligo[acc]=Hash.new
	# p d1.to_a.each_slice(7).to_a
	arrd1=d1.to_a.each_slice(7).to_a
	# p arrd1
	ptlen=arrd1.length
	arrd1.each_with_index{|part,ptidx|
		serlen=part.length
		easet=ptidx+1
		if ptidx == 0 then
			if part.length > 2 then
				part.each_with_index{|ol,idx|
					ono=ol[0].to_i
					seq=ol[1]
					unless d7aeoligo[acc][easet] then
						d7aeoligo[acc][easet]=Hash.new
					end
					if idx==0 then
						d7aeoligo[acc][easet]["#{ono}F"]=seq
					elsif idx==serlen-1 then
						if arrd1[ptidx+1] then
							if arrd1[ptidx+1].length > 2 then
								unless d7aeoligo[acc][easet+1] then
									d7aeoligo[acc][easet+1]=Hash.new
								end
								d7aeoligo[acc][easet+1]["#{ono}F"]=seq
								cmplseq=Bio::Sequence::NA.new(seq).reverse_complement.upcase
								d7aeoligo[acc][easet]["#{ono}R"]=cmplseq
							else
								cmplseq=Bio::Sequence::NA.new(seq).reverse_complement.upcase
								d7aeoligo[acc][easet]["#{ono}R"]=cmplseq
							end
						end
					else
						cmplseq=Bio::Sequence::NA.new(seq).reverse_complement.upcase
						d7aeoligo[acc][easet]["#{ono}R"]=cmplseq
					end
				}
			else	
				part.each_with_index{|ol,idx|
					ono=ol[0].to_i
					seq=ol[1]			
					cmplseq=Bio::Sequence::NA.new(seq).reverse_complement.upcase
					d7aeoligo[acc][easet]["#{ono}R"]=cmplseq
				}
			end
		elsif ptidx < ptlen-1 then
			if part.length > 2 then
				part.each_with_index{|ol,idx|
					ono=ol[0].to_i
					seq=ol[1]
					unless d7aeoligo[acc][easet] then
						d7aeoligo[acc][easet]=Hash.new
					end
					if idx==serlen-1 then
						if arrd1[ptidx+1] then
							if arrd1[ptidx+1].length > 2 then
								unless d7aeoligo[acc][easet+1] then
									d7aeoligo[acc][easet+1]=Hash.new
								end
								d7aeoligo[acc][easet+1]["#{ono}F"]=seq
								cmplseq=Bio::Sequence::NA.new(seq).reverse_complement.upcase
								d7aeoligo[acc][easet]["#{ono}R"]=cmplseq
							else
								cmplseq=Bio::Sequence::NA.new(seq).reverse_complement.upcase
								d7aeoligo[acc][easet]["#{ono}R"]=cmplseq
							end
						end
					else
						cmplseq=Bio::Sequence::NA.new(seq).reverse_complement.upcase
						d7aeoligo[acc][easet]["#{ono}R"]=cmplseq
					end
				}
			else	
				part.each_with_index{|ol,idx|
					ono=ol[0].to_i
					seq=ol[1]			
					cmplseq=Bio::Sequence::NA.new(seq).reverse_complement.upcase
					d7aeoligo[acc][easet]["#{ono}R"]=cmplseq
				}
			end
		else
			if part.length > 2 then
				unless d7aeoligo[acc][easet] then
					d7aeoligo[acc][easet]=Hash.new
				end
				part.each_with_index{|ol,idx|
					ono=ol[0].to_i
					seq=ol[1]

					cmplseq=Bio::Sequence::NA.new(seq).reverse_complement.upcase
					d7aeoligo[acc][easet]["#{ono}R"]=cmplseq

				}
			else	
				part.each_with_index{|ol,idx|
					ono=ol[0].to_i
					seq=ol[1]			
					cmplseq=Bio::Sequence::NA.new(seq).reverse_complement.upcase
					d7aeoligo[acc][easet-1]["#{ono}R"]=cmplseq
				}
			end
		end
	}
}
p d7aeoligo
CSV.open("aesoe2_oligo.csv","w") do |csv|
	d7aeoligo.each{|acc,d1|
		d1.each{|set,d2|
			d2.each{|i,seq|
				p "#{acc} #{set} #{i} #{seq}"
				csv << ["ae2_#{acc}_#{set}_#{i}", seq]
			}
		}
	}
end

File.open("aesoe2_oligo.fasta","w") do |file|
	d7aeoligo.each{|acc,d1|
		d1.each{|set,d2|
			d2.each{|i,seq|
				# p "#{acc} #{set} #{i} #{seq}"
				file.puts("ae2_#{acc}_#{set}_#{i}")
				file.puts(seq)
			}
		}
	}
end
