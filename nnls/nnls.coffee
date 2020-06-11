# Implements NNLS Algorithm as Described
# by RASMUS BRO and SIJMEN DE JONG
# Jorunal of Chemometrics, Vol 11 393-401 1997
# (c) by Jeremy Deuel, 2017
# jeremy@deuel.ch
# Uses numerical.js (www.numericaljs.com)


#Overload Array type
#transpose
Array.prototype.t = () ->
	if @.dim() == 2
		r = xoes([@.size(1),@.size(0)],NaN)
		for m in [0...@.size(0)]
			for n in [0...@.size(1)]
				r[n][m] = @[m][n]
		return r
	else if @.dim() == 1
		return @.map (x) -> [x]
# get size vector
Array.prototype.size = (a=-1) ->
	if a==-1
		numeric.dim(@)
	else
		numeric.dim(@)[a]
# get number of dimensions (scalar)
Array.prototype.dim = () ->
	numeric.dim(@).length
# get dot product
Array.prototype.dot = (b) ->
	numeric.dot(@,b)
# get max value
Array.prototype.max = () ->
	if (@.dim()>1)
		t = @.map (x) ->
			x.max()
		return t.max()
	m = @[0]
	for i in [0...@.length]
		m = @[i] if @[i]>m
	return m
# get min value
Array.prototype.min = () ->
	if (@.dim()>1)
		t = @.map (x) ->
			x.max()
		return t.max()
	m = @[0]
	for i in [0...@.length]
		m = @[i] if @[i]<m
	return m
# get index of max value
Array.prototype.argmax = () ->
	throw "can only calculate argmax for 1D-arrays" if @.dim() != 1
	m = @.max()
	for i in [0...@.length]
		return i if @[i] == m
#select values with binary vectors
Array.prototype.select = (m,n) ->
	if m instanceof Array
        if n instanceof Array
            throw("dimensions do not match: input: "+@.size()[0]+","+@.size()[1]+", m-array:"+m.length+" n-array:"+n.length) if @.size()[0] != m.length || @.size()[1] != n.length
            return @.select(m).t().select(n).t()
        else
            throw("dimensions do not match: input: "+@.size()[0]+", m-array:"+m.length) if @.size()[0] != m.length
            r = []
            for i in [0...m.length]
                r.push(@[i]) if m[i] 
            return r
    else
        if n instanceof Array
            throw("dimensions do not match: input: "+@.size()[1]+" n-array: "+n.length)if @.size()[1] != n.length
            return @.t().select(n).t()
        else
            return []
#set values with binary vectors
Array.prototype.set = (x,m,n) ->
	if m instanceof Array
        if n instanceof Array
            throw("dimensions do not match: target: "+@.size()[0]+","+@.size()[1]+" m-vector: "+m.length+" n-vector: "+n.length) if @.size()[0] != m.length || @.size()[1] != n.length
            throw("dimensions do not match: source: "+x.size()[0]+","+x.size()[1]+" m-vector: "+m.length+" n-vector: "+n.length) if x.size()[0] != m.length || x.size()[1] != n.length
            for i in [0...m.length]
                for j in [0...m.length]
                    @[i][j] = x[i][j] if m[i] && n[j]
            return @
        else
            throw("dimensions do not match: target: "+@.size()[0]+" m-vector: "+m.length) if @.size()[0] != m.length
            throw("dimensions do not match: source: "+x.size()[0]+" m-vector: "+m.length) if x.size()[0] != m.length
            for i in [0...m.length]
                @[i] = x[i] if m[i]
            return @
    else
        if n instanceof Array
            throw("dimensions do not match: target: "+@.size()[1]+" n-vector: "+n.length)if @.size()[1] != n.length
            throw("dimensions do not match: source: "+x.size()[1]+" n-vector: "+n.length)if x.size()[1] != n.length
            @.set(@.t().set(x.t(),n).t())
            return @
        else
        	#throw("dimensions do not match: target: "+@.size()[0]+","+@.size()[1]+", source: "+x.size()[0]+","+x.size()[1]) if x.size()[0] != @.size()[0] or x.size()[1] != @.size()[1]
            for i in [0...@.size()[0]]
                @[i] = x[i]
            return @
# negate
Array.prototype.not = () ->
    throw "can only bool inverse 1D arrays" if @.dim() != 1
    @.map (x) -> !x
# check if both arrays are 1D and of same length
Array.prototype.check1D = (x) ->
	throw "can only use 1D arrays" if @.dim() != 1 or x.dim() != 1
	throw "arrays dont have same length" if @.length != x.length
# rep y for dimensions defined in x
@xoes = (x,y) ->
	x = [x] if not (x instanceof Array)
	r = []
	if x.length>1
		for i in [0...x[0]]
			r.push(xoes(x.slice(1),y))
		return r
	for i in [0...x[0]]
		r.push(y)
	return r
# generate array with 0 of dimensions x
@zeros = (x) -> xoes(x,0)
# dito but with 1
@ones = (x) -> xoes(x,1)
# dito but with false
@falses = (x) -> xoes(x,false)
#dito but with true
@trues = (x) -> xoes(x,true)
#perform dotwise function
Array.prototype.check = (x,fun) ->
	x = xoes(@.size(),x) if not (x instanceof Array)
	@.check1D(x)
	@.map fun, x
Array.prototype.and = (x) -> @.check(x,(cV,i)->cV and @[i])
Array.prototype.or  = (x) -> @.check(x,(cV,i)->cV or @[i])
Array.prototype.lt  = (x) -> @.check(x,(cV,i)->cV < @[i])
Array.prototype.lteq= (x) -> @.check(x,(cV,i)->cV <= @[i])
Array.prototype.gt  = (x) -> @.check(x,(cV,i)->cV > @[i])
Array.prototype.gteq= (x) -> @.check(x,(cV,i)->cV >= @[i])
Array.prototype.add = (x) -> @.check(x,(cV,i)-> cV+@[i])
Array.prototype.sub = (x) -> @.check(x,(cV,i)-> cV-@[i])
Array.prototype.mul = (x) -> @.check(x,(cV,i)-> cV*@[i])
Array.prototype.div = (x) -> @.check(x,(cV,i)-> cV/@[i])
Array.prototype.pow = (x) -> @.check(x,(cV,i)-> Math.pow(cV,@[i]))
Array.prototype.mod = (x) -> @.check(x,(cV,i)-> cV%@[i])
#calculate pseudoinverse
Array.prototype.pinv = () ->
	z = numeric.svd(@)
	U = z.U
	S = z.S
	V = z.V
	m = @.length
	n = @[0].length
	tol = Math.max(m,n)*numeric.epsilon*z.S[0]
	M = S.length
	Sinv = new Array(M)
	for i in [M-1..0]
		if S[i]>tol
			Sinv[i] = 1/S[i]
		else
			Sinv[i] = 0
	return U.dot(numeric.diag(Sinv)).dot(V.t())
#sum up
Array.prototype.sum = () ->
	if @.dim()>1
		return @.map((x)->x.sum()).sum()
	else
		s = 0
		@.map((x)->s+=x)
		return s
# Returns the index of the i'th element, that is true
Array.prototype.ithTrue = (i) ->
	throw "ithTrue only works with one dimension" if @.dim()>1
	j = 0
	for k in [0...@.length]
		if @[k]
			if j==i
				return k
			else
				j+=1
	return null
#performs non-negative least squares on matrix Z and vector x, returns result vector d
@nnls =(Z, x,tol=1e-8) ->
	n = Z.size(0)
	P = falses(n)
	d = zeros(n)
	w = Z.dot(x.sub(Z.t().dot(d)))
	maxiter = 30*n
	itercount = 0
	while P.sum()<n and w.select(P.not()).gt(tol).sum()
		itercount += 1
		m = P.not().ithTrue(w.select(P.not()).argmax())
		P[m] = true
		s = zeros(n)
		ZP = Z.select(P)
		ZPT = Z.select(P).t()
		sp = ZP.dot(ZPT).pinv().dot(ZP).dot(x)
		j = 0
		for i in [0...P.length]
			if P[i]
				s[i] = sp[j]
				j+=1
		while sp.min()<0
			alpha = -1*(d.select(P).div(sp.sub(d.select(P))).min())
			d = d.add(s.sub(d).mul(alpha))
			P = P.and(d.gt(tol))
			s = zeros(n)
			ZP = Z.select(P)
			ZPT = Z.select(P).t()
			sp = ZP.dot(ZPT).pinv().dot(ZP).dot(x)
			j = 0
			for i in [0...P.length]
				if P[i]
					s[i] = sp[j]
					j+=1
		d = s
		w = Z.dot(x.sub(Z.t().dot(d)))
		if s[m] < tol
			break
		if itercount > maxiter
			break
	return d
#calc residue for nnls
@calcSpectrum = (Z,d) ->
	Z.t().dot(d)
@residueVector = (Z,x,d) ->
	x.sub(calcSpectrum(Z,d))
#calc residue (sqrt of sum of squares)
@residue = (Z,x,d) ->
	Math.sqrt(residueVector(Z,x,d).pow(2).sum())


class @Spectrum
	@load: (url) ->
		a = new Spectrum()
		a.loaded = false
		$.get url, null, (response, status, jqXHR) ->
			if status == 'success'
				a.parseCSV(response)
			else
				throw "Could not load file "+url
		return a
	@names: {}
	@selected: {}
	@loadNames: (callback) ->
		if typeof(Storage) != "undefined" then selection = JSON.parse(localStorage.getItem('references'))
		console.log(selection)
		if typeof(selection) != "object" or selection == null then selection = {}
		if Object.keys(@names).length == 0
			$.getJSON {url: "refspectra.json", cache: false}, (response, status, jqXHR) ->
				if status == 'success'
					#hasLocalStorage = (typeof(Storage) != "undefined")
					#if hasLocalStorage	
					#	v = localStorage.getItem("referenceVersion")
					#	if (v and v== response.version)
					#		Spectrum.names = JSON.parse(localStorage.getItem("referenceSpectra"))
					#		console.log("retrieving Spectra data, version "+v+" from local browser storage")
					#		console.log(Spectrum.names)
					#		callback(Spectrum.names)
					#		return
					#	else 
					#		console.log("local version "+v+" not identical with server version "+response.version)
					Spectrum.names = response.referenceSpectra
					for i in Object.keys(Spectrum.names)
						Spectrum.names[i].spectrum = Spectrum.load(Spectrum.names[i].url)
						Spectrum.names[i].selected = selection[i] if typeof(selection[i]) == "boolean"
					Spectrum.names['Precipitate'] = {
						selected: if typeof(selection[i]) == "boolean" then selection[i] else true
						url: ""
						spectrum: Spectrum([[200..800],ones(601).mul(0.01)])
						color: "#999999"
						}
					#saving local spectra only of pronto
					#localStorage.setItem("referenceVersion",response.version) if hasLocalStorage
					callback(Spectrum.names)
				else
					console.log("ERROR: Could not load refspectra.json (Default reference spectrum settings")		
	constructor: (data) ->
		if data instanceof Array
			throw "invalid input data for spectrum: not a 2D array" if data.dim() != 2
			if data.length != 2
				data = data.t()
			throw "invalid input data for spectrum: not two columns given" if data.length != 2
			a = new Spectrum()
			a.loaded = false
			a.loadArrays(data[0],data[1])
			return a
		if typeof(data)=="string" or data instanceof String
			a = new Spectrum()
			a.loaded = false
			a.parseCSV(data)
			return a
	loadArrays: (wl, ext) ->
		@data = {}
		for i in [0...wl.length]
			@data[Math.round(wl[i])] = ext[i]
		@wlmin = Math.round(wl.min())
		@wlmax = Math.round(wl.max())
		@wlstep = (@wlmax-@wlmin)/(wl.length-1)
		@extmin = ext.min()
		@extmax = ext.max()
		if @extmax>2 then @extmax = 2
		@loaded = true
	makeCSV: ->
		csvContent = "data:txt/csv;charset=utf-8,Wavelength (nm),Extinction\n"
		$.each @data, (x,y) ->
			csvContent += x+","+y+"\n"
		csvContent
	parseCSV: (data) ->
		lines = data.split("\r\n")
		lines = data.split("\n") if (lines.length<2)
		lines = data.split("\r") if (lines.length<2)
		throw "invalid input data: no lines detected: "+data if lines.length<2
		wl = []
		ext = []
		conc = 1
		if lines.length>2
			res = lines[1].match /^conc:\s(\d+[\.,]?\d*)/i
			if res
				conc = parseFloat(res[1])
		for l in lines
			res = l.match /^\s*"?(-?\d+\.?\d*)"?\s*[,;]\s*"?(-?\d+\.?\d*)"?/
			if res
				wl.push(parseFloat(res[1]))
				ext.push(parseFloat(res[2]))
			else if res==null and wl.length
				break
		if wl.length == ext.length and wl.length>0
			if (conc != 1)
				ext = ext.div(conc)
		@loadArrays(wl, ext)
	plotData: (conc=1,wlmin,wlmax)->
		data = []
		for i in Object.keys(@data).sort()
			continue if i<wlmin
			continue if i>wlmax
			data.push {x: i, y: @data[i]*conc}
			lasti = i
		return data
	interpolateWL: (wl) ->
		wl = Math.round(wl)
		return @data[wl] if @data[wl] != null
		return NaN if wl<@wlmin
		return NaN if wl>@wlmax
		r = (wl-@wlmin)/@wlstep
		rf = Math.floor(r)*@wlstep + @wlmin
		rc = Math.ceil(r)*@wlstep + @wlmin
		r = r%1
		return @data[rf]*(1-r)+@data[rc]*r #linear interpolation
	wl: (wlmin,wlmax=null,wlstep=null) ->
		if wlmax != null
			wlstep = 1 if wlstep ==null
			min = wlmin/wlstep
			max = wlmax/wlstep
			r = []
			for wl in [min..max]
				r.push wl*wlstep
			return @wl(r)
		if wlmin instanceof Array
			r = []
			for wl in wlmin
				r.push @interpolateWL(wl)
			return r
		if typeof(wlmin) == "number"
			return @interpolateWL(wlmin)
	dcnv: (min=null, max=null, step=null) ->
		refs = []
		names = []
		a = @dcnv
		$.each Spectrum.names, (name, value) ->
			if value.selBox.is(':checked')
				names.push name
				refs.push value.spectrum
		if refs.length < 1
			console.log("no dncv with less than two spectra selected")
			return null
		if min==null
			min = $.map refs, (x) -> x.wlmin
			min = [min.max(),@wlmin].max()
			min = [Spectrum.constraints.min.constraint,min].max() if Spectrum.constraints.min.constraint>0
		if max==null
			max = $.map refs, (x) -> x.wlmax
			max = [max.min(),@wlmax].min()
			max = [Spectrum.constraints.max.constraint,max].min() if Spectrum.constraints.max.constraint>0
		if step==null
			step = $.map refs, (x) -> x.wlstep
			step = [step.max(),@wlstep].max()
			step = [Spectrum.constraints.step.constraint,step].max() if Spectrum.constraints.step.constraint>0
		if (min>=max)
			console.log("no dncv with min="+min+">max="+max)
			return null
		ref_matrix = $.map refs, (x) -> [x.wl(min,max,step)]
		x_vector = @wl(min,max,step)
		sol = nnls(ref_matrix,x_vector)
		res = {'result':{},'residue':null,'residueVector':[],'params':{'min':min,'max':max,'step':step}}
		for i in [0...names.length]
			res['result'][names[i]] = sol[i]
		res['residue'] = residue(ref_matrix,x_vector,sol)
		cS = calcSpectrum(ref_matrix,sol)
		res['extmax'] = [cS.max(),x_vector.max()].max()
		rV = residueVector(ref_matrix,x_vector,sol)
		res['calcSpectrum'] = []
		res['calcSpectrum'].push {x: i*step, y: cS.pop()} for i in [max/step..min/step]
		res['residueVector'] = []
		res['residueVector'].push {x: i*step, y: rV.pop()} for i in [max/step..min/step]
		return res
currentSpectrum = null
plotConfigObject = (wlmin,wlmax) ->
	{
		data: 
			datasets: [{
				label:'Sample'
				pointRadius: 0
				data:currentSpectrum.plotData(1,wlmin,wlmax)
				fill: false
				backgroundColor: 'black'
				borderColor: 'black'
				borderWidth: 2
				}]
		options:
			title:
				display: true
				text: 'Deconvoluted Spectrum'
			scales:
				xAxes: [{
					ticks:
						min: if Spectrum.constraints.min.constraint then Spectrum.constraints.min.constraint else currentSpectrum.wlmin
						max: if Spectrum.constraints.max.constraint then Spectrum.constraints.max.constraint else currentSpectrum.wlmax
				}]
				yAxes: [{
					ticks:
						min: 0
						max: currentSpectrum.extmax*1.1
				}]
		}		
		
updateDcnv = () ->
	return if not (currentSpectrum instanceof Spectrum)
	for i in Object.keys(Spectrum.names)
		if not Spectrum.names[i].spectrum.loaded
			setTimeout(f, 100);
			console.log("not all references loaded, trying again in 100ms")
			return
	if not currentSpectrum.loaded
		setTimeout(f,100)
		console.log("currentSpectrum not loaded, trying again in 100ms")
		return
	cC = $('#canvas_container')
	cC.empty()
	c = $ '<canvas>'
	cC.append c
	c[0].style.width = '100%'
	c[0].style.height = '100%'
	c[0].width = c.offsetWidth
	c[0].height = c.offsetHeight
	res = currentSpectrum.dcnv()	
	if res
		d = plotConfigObject(res.params.min,res.params.max)
	else
		d = plotConfigObject(currentSpectrum.wlmin,currentSpectrum.wlmax)	
	for i in Object.keys(Spectrum.names)
		if res then v = res.result[i]
		if (typeof(v) == "number")
			Spectrum.names[i].valueCell.html v.toFixed(4)
			d.data.datasets.push {
				label: i+" ("+v.toFixed(2)+"µM)"
				pointRadius: 0
				data: Spectrum.names[i].spectrum.plotData(v,res.params.min,res.params.max)
				fill: false
				backgroundColor: Spectrum.names[i].color
				borderColor: Spectrum.names[i].color
				borderWidth: 1
			}
		
		else
			Spectrum.names[i].valueCell.html ''
	if res.extmax
		d.options.scales.yAxes[0].ticks.max = Math.ceil(res.extmax*5.1)/5
	if (res && res.calcSpectrum) then d.data.datasets.push {
			label: "Calculated Spectrum"
			data: res.calcSpectrum
			pointRadius: 0
			fill: false
			backgroundColor: "#ff0000"
			borderColor: "#ff0000"
			borderWidth: 2
		}
	if res then Spectrum.constraints.min.text.val res.params.min
	if res then Spectrum.constraints.max.text.val res.params.max
	if res then Spectrum.constraints.step.text.val res.params.step
	if res && typeof(res.residue) == "number"
		Spectrum.residue.html (res.residue).toFixed(3)
	else
		Spectrum.residue.html ''
	c = Chart.Scatter(c[0].getContext("2d"),d)
updateConstraints = (evt) ->
	persistentConstraint = {}
	for i in Object.keys(Spectrum.constraints)
		v = parseFloat(Spectrum.constraints[i].text.val())
		if Spectrum.constraints[i].check.is(':checked')
			Spectrum.constraints[i].constraint = v
		else
			Spectrum.constraints[i].constraint = NaN
		if isNaN(v)
			Spectrum.constraints[i].text.val 'auto'
		else
			Spectrum.constraints[i].text.val v.toFixed(0)
		persistentConstraint[i] = {constraint: Spectrum.constraints[i].constraint}
	localStorage.setItem('constraints',JSON.stringify(persistentConstraint)) if typeof(Storage) != "undefined"
	updateDcnv()		
updateReferences = ->
	updateDcnv()
	if typeof(Storage) != "undefined"
		selected = {}
		for i in Object.keys(Spectrum.names)
			selected[i] = Spectrum.names[i].selBox.is(':checked')
			if selected[i]
				Spectrum.names[i].selBox.parent().parent().removeClass('deselectedRow')
			else
				Spectrum.names[i].selBox.parent().parent().addClass('deselectedRow')

		localStorage.setItem('references',JSON.stringify(selected))
buildReferences = ->
	rB = $('#referenceBox')
	rB.empty()
	table = $ '<table class="table table-sm">'
	table.append ($ '<thead><tr><th class="firstCol">&nbsp;</th><th>Substance</th><th>Conc (µM)</th></tr></thead>')
	rB.append table
	tbody = $ '<tbody>'
	table.append tbody 
	for i in Object.keys(Spectrum.names)
		if not Spectrum.names[i].spectrum.loaded
			rB.append '<div class="alert alert-danger">Waiting for reference spectra to load...</div>'
			setTimeout(buildReferences, 100);
			console.log "Spectrum "+i+" not yet loaded, trying again in 100ms"
			return
		tr = $ '<tr>'
		tbody.append tr
		td = $ '<td class="firstCol">'
		tr.append td
		if Spectrum.names[i].selected
			cb = $ "<input type=\"checkbox\" checked=\"checked\"/>"
		else
			cb = $ "<input type=\"checkbox\" />" 
			tr.addClass('deselectedRow')
		Spectrum.names[i].selBox = cb
		td.append cb
		cb.data('ref',i)
		cb.change(updateReferences)
		td = $ '<td>'
		tr.append td
		td.html i
		td = $ '<td>'
		tr.append td
		Spectrum.names[i].valueCell = td
	Spectrum.constraints = {min: {}, max: {}, step: {}} if typeof(Spectrum.constraints) != "object" or Spectrum.constraints==null
	for i in Object.keys(Spectrum.constraints)
		tr = $ '<tr>'
		tbody.append tr
		td = $ '<td class="firstCol">'
		cb = $ "<input type=\"checkbox\">"
		Spectrum.constraints[i].check = cb
		cb.attr('checked','checked') if (Spectrum.constraints[i].constraint)>0
		cb.change updateConstraints
		td.append cb
		tr.append td
		tr.append $ '<td>wavelength '+i+'</td>'
		td = $ '<td>'
		tr.append td
		cb = $ "<input type=\"text\" value=\"auto\" placeholder=\"auto\" size=\"6\">"
		cb.val Spectrum.constraints[i].constraint.toFixed(0) if (Spectrum.constraints[i].constraint)>0
		Spectrum.constraints[i].text = cb
		cb.change updateConstraints
		td.append cb
	# residue
	Spectrum.residue = $ '<td>'
	tr = $ '<tr>'
	tbody.append tr
	tr.append '<td class="firstCol">&nbsp;</td>'
	tr.append '<td>Residue</td>'
	tr.append Spectrum.residue
	updateDcnv()
$ ->
	Spectrum.constraints = JSON.parse(localStorage.getItem('constraints')) if typeof(Storage) != "undefined"
	Spectrum.loadNames buildReferences
	
	if (location.hash && location.hash.length>10)
		try
			currentSpectrum = new Spectrum()
			currentSpectrum.parseCSV(atob(location.hash.substr(1)))
			throw "no valid string" if not currentSpectrum.loaded
			console.log("loaded test spectrum from url")
		catch e
			currentSpectrum = null
			console.log("failed to read a base64 encoded string: "+e)
	fdz = document.getElementById("file_drop_zone")
	fdz.addEventListener 'dragover',((evt) -> 
		evt.stopPropagation()
		evt.preventDefault()
		evt.dataTransfer.dropEffect = 'copy'
		evt.dataTransfer.effectAllowed = 'copy'
		),false
			
	fdz.addEventListener 'drop', ((evt) ->
		evt.stopPropagation()
		evt.preventDefault()
		dt = evt.dataTransfer
		if (dt.files && dt.files.length >0)
			fi = dt.files[0]
			reader = new FileReader()
			reader.onload = (evt) ->
				try
					currentSpectrum = new Spectrum()
					currentSpectrum.parseCSV reader.result
					updateDcnv()
				catch e
					alert "Failed to parse CSV for dropped file "+fi.name+": "+e
			reader.readAsText fi
			return
		),false
	fdz.addEventListener 'dragend', ((evt) ->
		dt = evt.dataTransfer
		if (dt.items)
			for i in [0...dt.items.length]
				dt.items.remove(i)
		else
			evt.dataTransfer.clearData()
		),false
	fdz.addEventListener 'click', ((evt) ->
		return if not (currentSpectrum instanceof Spectrum)
		link = document.createElement 'a'
		link.download = "spectrum.csv"
		link.href = encodeURI(currentSpectrum.makeCSV())
		link.click()
		), false