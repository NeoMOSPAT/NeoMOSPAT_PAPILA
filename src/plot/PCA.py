# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib import gridspec, colors as colo
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans, MeanShift, estimate_bandwidth
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler 

import src.IncludeFile as IncF

###############################################
## Plot All things related PCA
def plot_PCA(dff, variables, filename, showplot=False):
	df=dff.copy()
	df=df.sort_values(by='Deviation')
	df["Hora"]=df.index.hour

	###########
	## Do PCA
	variances, components, df = getPCA(df, variables)
	

	## If there are at least two variables 
	if len(variables)>2:
		###########
		## Plots PCA
		## Definitions for the subplots
		left, width = 0.07, 0.23
		bottom, height = 0.145, .72
		left_2 = left+width+0.07
		left_3 = left_2+width+0.01
		##  Adjustments so that bottom aligns
		adj1=0.234
		adj3=-0.09
		if len(variables)==4:
			adj1=adj1
		elif len(variables)==3:
			adj1=1.04*adj1
		elif len(variables)==5:
			adj1=1.02*adj1
		elif len(variables)==6:
			adj1=1.035*adj1
		elif len(variables)==7:
			adj1=1.042*adj1
		## Subploy details
		b_1 = [left, bottom+adj1, width, height-adj1]
		b_2 = [left_2, bottom, width, height]
		b_3 = [left_3, bottom+adj3, 0.38, height-adj3]

		################
		## Figure hh
		fig = plt.figure(figsize=(15,4.9))
		ax1 = fig.add_axes(b_1)
		ax2 = fig.add_axes(b_2)
		ax3 = fig.add_axes(b_3, projection='3d')
		##################
		BarChart(variables, variances, components, filename+"__PCA-BarChart", ax=ax1, showplot=showplot)
		plot_2D_PCA_hora(df, filename+"__PCA-2D", ax=ax2, showplot=showplot)
		plot_3D_PCA_hora(df, filename+"__PCA-3D", ax=ax3, showplot=showplot)
		plt.savefig(IncF.savepath+filename+'__PCA__hh.pdf', format="pdf")
		plt.close(fig)

		##################
		##################
		## Figure sd
		fig = plt.figure(figsize=(15,4.9))
		ax1 = fig.add_axes(b_1)
		ax2 = fig.add_axes(b_2)
		ax3 = fig.add_axes(b_3, projection='3d')
		##################
		BarChart(variables, variances, components, filename+"__PCA-BarChart", ax=ax1, showplot=showplot)
		plot_2D_PCA_sd(df, filename+"__PCA-2D", ax=ax2, showplot=showplot)
		plot_3D_PCA_sd(df, filename+"__PCA-3D", ax=ax3, showplot=showplot)
		plt.savefig(IncF.savepath+filename+'__PCA__sd.pdf', format="pdf")
		plt.close(fig)

		##################
		plot_3D_PCA_sd(df, filename+"__PCA-3D", ax=None,  showplot=showplot)
		plot_3D_PCA_hora(df, filename+"__PCA-3D", ax=None, showplot=showplot)

	##################
	plot_2D_PCA_sd(df, filename+"__PCA-2D", ax=None, showplot=showplot)
	plot_2D_PCA_hora(df, filename+"__PCA-2D", ax=None, showplot=showplot)
	BarChart(variables, variances, components, filename+"__PCA-BarChart", ax=None, showplot= showplot) 

	return None


########################################################
########################################################
########################################################
### Obtain PCA analysis
def getPCA(dff, features):
	df=dff.copy()
	## Separating out the variables to study and standaring the data
	df=df.dropna(how="any")
	x=df.loc[:, features].values
	dates=df.index
	df.reset_index(drop=True, inplace=True)
	x=StandardScaler().fit_transform(x)	

	## Get number of variables, hence PC what needs to be done
	num_dim=len(features)

	## PCA
	pca=PCA(n_components=num_dim)
	principalComponents=pca.fit_transform(x)
	variances=pca.explained_variance_ratio_
	components=[[ '%.2f' % elem for elem in my_list ] for my_list in pca.components_]

	## Data frame with princiapl components
	cols=["PC%s"%(i+1) for i in range(num_dim)]
	principalDf=pd.DataFrame(data = principalComponents, columns = cols)

	finalDf=pd.concat([principalDf, df], axis = 1)
	finalDf=finalDf.loc[:,~finalDf.columns.duplicated()]
	finalDf.index=dates
	return [variances, components, finalDf]


##############################################
## Bar Chart with the % of information per Principal components and how it's composed of
def BarChart(variables, variances, components, filename, ax=None, showplot=False):

	## Number of variables/PC
	num_dim=len(variances)
	## Get the strings for the y axis [PC1, PC2, ....]
	objects = ["PC"+str(i+1) for i in range(num_dim)]
	y_pos = np.arange(len(objects))


	#########
	## Figure
	#fig = plt.figure(figsize=(4,4))
	#ax = fig.add_subplot(111)	
	## Check if new or old plot
	BoolOwn=False
	if ax==None:
		BoolOwn=True
	ax=ax or plt.gca()
	## Plot bar plot
	ax.bar(y_pos, variances*100, align='center', alpha=0.5)

	ax.set_title('Principal Components Variance')	
	ax.set_ylabel("%")
	ax.set_yticks(np.arange(0, 100, 20))
	ax.set_xticks(y_pos)
	ax.set_xticklabels(objects)

	ax.yaxis.grid(linestyle='--') 
	## Change x axis given how many dimensions and put 0-100 %
	ax.axis([-0.6, num_dim-0.4, 0, 100])

	#########
	## Table
	font_size=8
	bbox=[0, -0.67, 1, 0.47]
	mpl_tabl = ax.table(cellText = components, rowLabels = variables, bbox=bbox, rowLoc="left", cellLoc="center", colLabels=objects, loc='center', rowColours=['lightblue']*num_dim, colColours=['lightblue']*num_dim)
	mpl_tabl.auto_set_column_width(col=num_dim)
	mpl_tabl.auto_set_font_size(False)
	mpl_tabl.set_fontsize(font_size)
	plt.subplots_adjust( bottom=0.4, left=0.3)#left=0.3

	## Save File
	plt.savefig(IncF.savepath+filename+'.pdf', format="pdf")
	#fig.savefig(savepath+filename+'.pdf', format="pdf")
	
	
	## Show plot
	if BoolOwn:
		if showplot :
			plt.show()
			return None
		else:
			plt.clf()
			plt.close()
			#plt.close(fig)
			return None
	else:
		return None


###############################################
## Plot PCA with hour being heighlighted
def plot_3D_PCA_hora(df, filename, ax=None, showplot=False):
	#df=df.sort_values(by='PMG')
	#####
	## Figure
	#fig = plt.figure()
	#ax = fig.add_subplot(111, projection='3d')
	BoolOwn=False
	if ax==None:
		BoolOwn=True
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

	################
	# Marker Size
	size=18
	## Marker Color
	cmap = plt.get_cmap('hsv', 24)
	## Edge Colors
	my_cmap=np.array([cmap(i) for i in np.linspace(0, 1, 24)])
	my_cmap[:,0:3] *= 0.8

	#################
	## Obtain edge colors
	lista=df["Hora"].values 
	lista2=np.zeros(shape=(len(lista),4))
	key = [[i, my_cmap[i]]for i in range(24)]
	for k in key:
		lista2[lista == k[0]] = k[1]
	mcolors=lista2
 
 	#############
	## Plot
	p=ax.scatter(df["PC1"], df["PC2"], df["PC3"], s=size, cmap=cmap, c=df["Hora"], edgecolors=mcolors, linewidth=0.5, alpha=1)
	ax.set_xlabel("PC1")
	ax.set_ylabel("PC2")
	ax.set_zlabel("  PC3")

	#if BoolOwn==True:

	## Add colorbar
	cbar = plt.colorbar(p, pad=0.12)
	## Specify which points will have a legend
	cbar.set_ticks(0.45+np.arange(24)*24.0/25)
	cbar.set_ticklabels(["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23"])
	#cbar_ax = fig.axes[-1]
	#cbar_ax.tick_params(labelsize=9)

	ax.set_title('Three Component PCA')	
	azim=-56.82459677419354
	elev=21.720779220779264 
	ax.view_init(elev, azim)
	#plt.tight_layout()

	#plt.subplots_adjust(right=0.9)#left=0.3
	
	if BoolOwn:
		## Save Figure
		plt.savefig(IncF.savepath+filename+'__hh.pdf', format="pdf",bbox_inches='tight')
		## Show plot
		if showplot:
			plt.show()
			#print('ax.azim {}'.format(ax.azim))
			#print('ax.elev {}'.format(ax.elev))
			return None
		else:
			plt.clf()
			plt.close()
			return None
	else:
		return None


###############################################
## Plot PCA with Deviation being heighlighted
def plot_3D_PCA_sd(df, filename, ax=None, showplot=False):
	df=df.sort_values(by='PMG')

	#####
	## Figure
	#fig = plt.figure()
	#ax = fig.add_subplot(111, projection='3d')
	BoolOwn=False
	if ax==None:
		BoolOwn=True
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

	# Marker Size
	size=18
	## Marker Color
	cmap = plt.get_cmap('YlGnBu', 5)
	## Edge Colors
	my_cmap=np.array([cmap(i) for i in np.linspace(0, 1, 5)])
	my_cmap[:,0:3] *= 0.8

	## Obtaining lists of colors
	colordata=df["Deviation"]
	num0=len(colordata[colordata==0])
	num1=len(colordata[colordata==1])
	num2=len(colordata[colordata==2])
	num3=len(colordata[colordata==3])
	num4=len(colordata[colordata==4])

	my_cmap[3]=np.array([17.0/255, 62.0/255, 105.0/255,1])
	my_cmap[4]=np.array([9.0/255, 19.0/255, 71.0/255,1])

	mcolors=[my_cmap[0]]*num0+[my_cmap[1]]*num1+[my_cmap[2]]*num2+[my_cmap[3]]*num3+[my_cmap[4]]*num4


	## Plot
	p=ax.scatter(df["PC1"], df["PC2"], df["PC3"], s=size, cmap=cmap, c=colordata, edgecolors=mcolors,  vmin=-1, vmax=4.5, linewidth=0.5, alpha=1)
	ax.set_xlabel("PC1")
	ax.set_ylabel("PC2")
	ax.set_zlabel("  PC3")

	## Add colorbar
	cbar = plt.colorbar(p, pad=0.12)
	## Specify which points will have a legend
	cbar.set_ticks([-0.45, 0.65, 1.75, 2.78, 3.95])
	## Change label
	cbar.set_ticklabels(['$\pm 0 \sigma$', '$\pm 1 \sigma$', '$\pm 2 \sigma$', '$\pm 3 \sigma$', '$\pm 4 \sigma$'])

	ax.set_title('Three Component PCA')	
	azim=-56.82459677419354
	elev=21.720779220779264
	ax.view_init(elev, azim)
	#fig.tight_layout()


	
	if BoolOwn:
		## Save Figure
		plt.savefig(IncF.savepath+filename+'__sd.pdf', format="pdf", bbox_inches='tight')
		## Show plot
		if showplot:
			return plt.show()
		else:
			plt.close()
			return None
	else:
		return None


###############################################
## Plot PCA with Deviation being heighlighted
def plot_2D_PCA_hora(df, filename, ax=None, showplot=False):
	df=df.sort_values(by='PMG')
	#####
	## Figure
	#fig = plt.figure()
	#ax = fig.add_subplot(111)
	BoolOwn=False
	if ax==None:
		BoolOwn=True
	ax=ax or plt.gca()

	################
	# Marker Size
	size=18
	## Marker Color
	cmap = plt.get_cmap('hsv', 24)
	## Edge Colors
	my_cmap=np.array([cmap(i) for i in np.linspace(0, 1, 24)])
	my_cmap[:,0:3] *= 0.8

	#################
	## Obtain edge colors
	lista=df["Hora"].values 
	lista2=np.zeros(shape=(len(lista),4))
	key=[[i, my_cmap[i]]for i in range(24)]
	for k in key:
		lista2[lista==k[0]]=k[1]
	mcolors=lista2

	###################
	## Plot
	p=ax.scatter(df["PC1"], df["PC2"], s=size, cmap=cmap, c=df["Hora"], edgecolors=mcolors, linewidth=0.5, alpha=1)	
	ax.set_xlabel("PC1")
	ax.set_ylabel("PC2")

	if BoolOwn==True:
		## Add colorbar
		cbar = plt.colorbar(p)
		## Specify which points will have a legend
		cbar.set_ticks(0.45+np.arange(24)*24.0/25)
		cbar.set_ticklabels(["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23"])
		#cbar_ax = fig.axes[-1]
		#cbar_ax = plt.axes[-1]
		#cbar_ax.tick_params(labelsize=9)

	ax.set_title('Two Component PCA')	
	#plt.tight_layout()
	#fig.tight_layout()



	if BoolOwn:
		## Save Figure
		#fig.savefig(savepath+filename+'__hh.pdf', format="pdf")
		plt.savefig(IncF.savepath+filename+'__hh.pdf', format="pdf",bbox_inches='tight')

		## Show plot
		if showplot:
			plt.show()
			return None
		else:
			#plt.close(fig)
			plt.clf()
			plt.close()
			return None
	else:
		return None


###############################################
## Plot PCA with Deviation being heighlighted
def plot_2D_PCA_sd(df, filename, ax=None, showplot=False):

	#####
	## Figure
	#fig = plt.figure()
	#ax = fig.add_subplot(111)
	BoolOwn=False
	if ax==None:
		BoolOwn=True
	ax=ax or plt.gca()


	# Marker Size
	size=18
	## Marker Color
	cmap = plt.get_cmap('YlGnBu', 5)
	## Edge Colors
	my_cmap=np.array([cmap(i) for i in np.linspace(0, 1, 5)])
	my_cmap[:,0:3] *= 0.8

	## Obtaining lists of colors
	colordata=df["Deviation"]
	num0=len(colordata[colordata==0])
	num1=len(colordata[colordata==1])
	num2=len(colordata[colordata==2])
	num3=len(colordata[colordata==3])
	num4=len(colordata[colordata==4])

	my_cmap[3]=np.array([17.0/255, 62.0/255, 105.0/255,1])
	my_cmap[4]=np.array([9.0/255, 19.0/255, 71.0/255,1])

	mcolors=[my_cmap[0]]*num0+[my_cmap[1]]*num1+[my_cmap[2]]*num2+[my_cmap[3]]*num3+[my_cmap[4]]*num4

	## Plot
	p=ax.scatter(df["PC1"], df["PC2"], s=size, cmap=cmap, c=colordata, edgecolors=mcolors,  vmin=-1, vmax=4.5, linewidth=0.5, alpha=1)
	ax.set_xlabel("PC1")
	ax.set_ylabel("PC2")

	if BoolOwn==True:
		## Add colorbar
		cbar = plt.colorbar(p)
		## Specify which points will have a legend
		cbar.set_ticks([-0.45, 0.65, 1.75, 2.78, 3.95])
		## Change label
		cbar.set_ticklabels(['$\pm 0 \sigma$', '$\pm 1 \sigma$', '$\pm 2 \sigma$', '$\pm 3 \sigma$', '$\pm 4 \sigma$'])

	ax.set_title('Two Component PCA')	
	#fig.tight_layout()

	if BoolOwn:
		## Save Figure
		plt.savefig(IncF.savepath+filename+'__sd.pdf', format="pdf",bbox_inches='tight')
		## Show plot
		if showplot:
			return plt.show()
		else:
			plt.close()
			return None
	else:
		return None


###############################################
## K-means clustering 2D
def PC1_PC2_kmeans(finalDf, filename, showplot):
	## K-means Clustering
	kmeans = KMeans(n_clusters=2)
	X=finalDf[["PC1", "PC2"]].values
	kmeans.fit(X)
	## List of numbers of what cluster they belong to for cmap in scatter plot
	y_kmeans = kmeans.predict(X)
	## Cluster centers
	centers = kmeans.cluster_centers_
		
	## Add column to which cluster it belongs to
	finalDf.loc[:,'Labels'] = kmeans.labels_
	finalDf=finalDf.sort_values(by='Labels')
	df0=finalDf.loc[finalDf['Labels']==0]
	df1=finalDf.loc[finalDf['Labels']==1]

	## Deleting Columns
	df0=df0[df0.columns.difference(["PC3", "PC4", "PC5","Ano","Mes", "Dia", "Labels"])]
	df1=df1[df1.columns.difference(["PC3", "PC4", "PC5","Ano","Mes", "Dia", "Labels"])]
	num_dim=len(df0.columns)

	## Getting Info
	mean_d0=pd.DataFrame([df0.mean()])
	mean_d1=pd.DataFrame([df1.mean()])
	sd_d0=pd.DataFrame([df0.std()])
	sd_d1=pd.DataFrame([df1.std()])
	min_d0=pd.DataFrame([df0.min()])
	min_d1=pd.DataFrame([df1.min()])
	max_d0=pd.DataFrame([df0.max()])
	max_d1=pd.DataFrame([df1.max()])
	rows=["Mean", "Standard Dev", "Min Value", "Max Value"]
	df_d0=pd.concat([mean_d0, sd_d0, min_d0, max_d0])
	df_d1=pd.concat([mean_d1, sd_d1, min_d1, max_d1])

	n0=len(df0)
	n1=len(df1)
			
	fig = plt.figure(figsize=(12,4))
	gs = gridspec.GridSpec(1, 2, width_ratios=[1, 2]) 
	
	###
	ax1 = plt.subplot(gs[0])
	ax1.scatter(df0["PC1"], df0["PC2"], color="orange", edgecolor="black")
	ax1.scatter(df1["PC1"], df1["PC2"], color="steelblue", edgecolor="black")

	###
	ax2 = plt.subplot(gs[1])
	vals0=df_d0.values
	vals0=[[ '%.2f' % elem for elem in my_list ] for my_list in vals0]
	vals1=df_d1.values
	vals1=[[ '%.2f' % elem for elem in my_list ] for my_list in vals1]

	columns=df0.columns

	sf=0.45

			
	bbox=[0.12, 0.55, 0.9, 0.45]
	mpl_tabl0 = ax2.table(cellText=vals0, colLabels=columns, rowLabels=rows, loc='center', rowColours=['orange']*len(rows), colColours=['orange']*num_dim, bbox=bbox)

	bbox=[0.12, 0.01, 0.9, 0.45]
	mpl_tabl1 = ax2.table(cellText=vals1, colLabels=columns, rowLabels=rows, loc='center', rowColours=['steelblue']*len(rows), colColours=['steelblue']*num_dim, bbox=bbox)
			

	mpl_tabl0.auto_set_font_size(False)
	mpl_tabl1.auto_set_font_size(False)
	font_size=9
	mpl_tabl0.set_fontsize(font_size)
	mpl_tabl1.set_fontsize(font_size)
	mpl_tabl0.auto_set_column_width(col=num_dim)
	mpl_tabl1.auto_set_column_width(col=num_dim)
	ax2.axis('off')	


	## Save Figure
	fig.savefig(IncF.savepath+filename+'__PCA-2__PC1-PC2__k-means.pdf', format="pdf")
	
	
	## Show plot
	if showplot:
		return plt.show()
	else:
		plt.close(fig)
		return None
