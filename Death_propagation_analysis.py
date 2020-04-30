import scipy.io
import numpy as np
import seaborn as sns
from scipy.spatial import Voronoi
import matplotlib
import matplotlib.pyplot as plt
import os
import math
import pandas as pd
from statistics import mode
import pysal as ps
import matplotlib as mpl
from scipy.stats.stats import pearsonr
from sklearn.neighbors import NearestNeighbors
import random


class Analysis:
    def __init__(self, file_path, output_path, time_frame, cell_line, treatment, signal_analyzed, n_scramble=1000, draw=False):
        self.file_path = '/'.join(file_path.split('/')[0:-1]) + '/'
        file_name = file_path.split('/')[-1]
        print(file_name)
        self.file_name = file_name.replace(".csv", "")
        data = pd.read_csv(file_path)
        new_path = output_path + "Results/" + self.file_name
        if not os.path.exists(new_path):
            os.makedirs(new_path)
        self.data_path = new_path + "/Data"
        self.pic_path = new_path + "/Pic"
        if not os.path.exists(self.data_path):
            os.makedirs(self.data_path)
        if not os.path.exists(self.pic_path):
            os.makedirs(self.pic_path)
        self.mean_time_death = None
        self.time_frame = time_frame
        self.cell_line = cell_line
        self.treatment = treatment
        self.signal_analyzed = signal_analyzed
        self.treatment_type = self.file_name[:self.file_name.find("_")]
        self.n_scramble = n_scramble
        self.listX = data["X"]
        self.lisY = data["Y"]
        self.n_instances = len(self.listX)
        self.die_times = data["times_of_death"]
        self.num_blobs = [x for x in range(max(self.die_times)+1)]
        self.experiment_time = len(self.num_blobs)
        self.XY = np.column_stack((self.listX, self.lisY))
        self.scrambles = self.create_scramble(self.XY, n=n_scramble)
        self.neighbors_difference_death_times = self.get_neighbors_difference_death_times()
        self.original_difference_death_times = self.neighbors_difference_death_times[0]
        self.scramble_signficance_95 = None
        self.scramble_signficance_98 = None
        self.scrample_mean_time_death = None
        self.statistic_score = self.assess_mean()
        self.death_perc_by_time = self.calc_death_perc()
        self.neighbors_stats = []
        self.difference_from_leader = []
        self.get_distance_time_from_leader()
        self.neighbors_list = []
        self.neighbors_list2 = []
        self.neighbors_list3 = []
        self.get_neighbors()
        self.who_is_leader = None
        self.death_wave = None
        self.identify_leader()
        self.draw() if draw else None
        self.dataframe = None
        self.create_dataframe()


    def create_dataframe(self):
        instances = []
        x_array = []
        y_array = []
        time_from_leader = []
        distance_from_leader = []
        first_neighbor_distance = []
        second_neighbor_distance = []
        third_neighbor_distance = []
        first_neighbor_time = []
        second_neighbor_time = []
        third_neighbor_time = []
        death_perc = []
        neighbors_number = []
        dead_neighbors_number = []
        has_1_dead_neighbors = []
        has_2_dead_neighbors = []
        has_3_dead_neighbors = []
        has_4_dead_neighbors = []
        has_5_dead_neighbors = []
        has_6_dead_neighbors = []
        has_7_dead_neighbors = []
        has_8_dead_neighbors = []
        has_9_dead_neighbors = []
        has_10_dead_neighbors = []
        cell_line = []
        treatment = []
        signal_analyzed = []
        file_name = []
        neighbors_distance_from_leader = []
        neighbors_stats=self.get_neighbors_stats(self.neighbors_list)
        second_neighbors_num = []
        avg_distance_from_neighbors = []
        for i in range(self.n_instances):
            instances.append(i)
            file_name.append(self.file_name)
            cell_line.append(self.cell_line)
            treatment.append(self.treatment)
            x_array.append(self.XY[i][0])
            y_array.append(self.XY[i][1])
            signal_analyzed.append(self.signal_analyzed)
            second_neighbors_num.append(len(self.neighbors_list2[i]))
            time_from_leader.append(self.difference_from_leader[i][0])
            distance_from_leader.append(self.difference_from_leader[i][1])
            neighbors_distance_from_leader.append(4)
            n = 0
            d = 0
            dis = 0
            for x in neighbors_stats[i]:
                n += 1
                dis += x[2]
                if x[1] >= 0:
                    d += 1
                    if d == 1:
                        first_neighbor_distance.append(x[2])
                        first_neighbor_time.append(x[1])
                    elif d == 2:
                        second_neighbor_distance.append(x[2])
                        second_neighbor_time.append(x[1])
                    elif d == 3:
                        third_neighbor_distance.append(x[2])
                        third_neighbor_time.append(x[1])
            n = 1 if n == 0 else n
            death_perc.append(d/n)
            avg_distance_from_neighbors.append(dis / n)
            neighbors_number.append(n)
            dead_neighbors_number.append(d)
            has_2_dead_neighbors.append(1 if d > 1 else 0)
            has_3_dead_neighbors.append(1 if d > 2 else 0)
            has_4_dead_neighbors.append(1 if d > 3 else 0)
            has_5_dead_neighbors.append(1 if d > 4 else 0)
            has_6_dead_neighbors.append(1 if d > 5 else 0)
            has_7_dead_neighbors.append(1 if d > 6 else 0)
            has_8_dead_neighbors.append(1 if d > 7 else 0)
            has_9_dead_neighbors.append(1 if d > 8 else 0)
            has_10_dead_neighbors.append(1 if d > 9 else 0)
            has_1_dead_neighbors.append(1 if d > 0 else 0)
            if d<  3:
                if d <= 2:
                    third_neighbor_distance.append(None)
                    third_neighbor_time.append(None)
                if d <= 1:
                    second_neighbor_distance.append(None)
                    second_neighbor_time.append(None)
                if d == 0:
                    first_neighbor_distance.append(None)
                    first_neighbor_time.append(None)
        for x in self.neighbors_list3[0]:
            neighbors_distance_from_leader[x] = 3
        for x in self.neighbors_list2[0]:
            neighbors_distance_from_leader[x] = 2
        for x in self.neighbors_list[0]:
            neighbors_distance_from_leader[x] = 1
        neighbors_distance_from_leader[0] = 0
        self.get_distance_from_local_leader(dead_neighbors_number)
        feature_table = pd.DataFrame(
            {'id': instances,
             'file_name': file_name,
             'cell_line': cell_line,
             'treatment': treatment,
             'signal_analyzed': signal_analyzed,
             'x_loc': x_array,
             'y_loc': y_array,
             'die_time': self.die_times,
             'die_time_real': [x * self.time_frame for x in self.die_times],
             'distance_from_leader':distance_from_leader,
             'time_from_leader':time_from_leader,
             'death_perc':death_perc,
             'first_neighbor_distance':first_neighbor_distance,
             'first_neighbor_time':first_neighbor_time,
             'second_neighbor_distance': second_neighbor_distance,
             'second_neighbor_time': second_neighbor_time,
             'third_neighbor_distance': third_neighbor_distance,
             'third_neighbor_time': third_neighbor_time,
             # 'total_density':self.density,
             # 'point_density':self.density_points,
             'neighbors_distance_from_leader':neighbors_distance_from_leader,
             'has_1_dead_neighbors': has_1_dead_neighbors,
             'has_2_dead_neighbors': has_2_dead_neighbors,
             'has_3_dead_neighbors': has_3_dead_neighbors,
             'has_4_dead_neighbors': has_4_dead_neighbors,
             'has_5_dead_neighbors': has_5_dead_neighbors,
             'has_6_dead_neighbors': has_6_dead_neighbors,
             'has_7_dead_neighbors': has_7_dead_neighbors,
             'has_8_dead_neighbors': has_8_dead_neighbors,
             'has_9_dead_neighbors': has_9_dead_neighbors,
             'has_10_dead_neighbors': has_10_dead_neighbors,
             'neighbors_number': neighbors_number,
             'second_neighbors_num': second_neighbors_num,
             'avg_distance_from_neighbors': avg_distance_from_neighbors,
             'dead_neighbors_number': dead_neighbors_number
             })

        location = self.data_path + "\File_{}.csv".format(self.file_name)
        feature_table.to_csv(location)
        self.dataframe = feature_table

    def create_scramble(self, xy, n=1000):
        scrambles = []
        for i in range(n):
            temp_copy = xy.copy()
            np.random.shuffle(temp_copy)
            scrambles.append(temp_copy)
        return scrambles

    def get_neighbors_difference_death_times(self):
        times = []
        times.append(self.get_time_from_neighbors(self.die_times, self.XY))
        for i in range(self.n_scramble):
            times.append(self.get_time_from_neighbors(self.die_times, self.scrambles[i]))
        return times

    def get_time_from_neighbors(self, times, points):
        vor = Voronoi(points)
        neighbors = vor.ridge_points
        t = []
        for i in range(len(neighbors)):
            t.append(abs(times[neighbors[i][0]] - times[neighbors[i][1]]))
        return t

    def get_neighbors(self):
        vor = Voronoi(self.XY)
        neighbors = vor.ridge_points
        leaders = []
        for i in range(self.n_instances):
            self.neighbors_list.append([])
            self.neighbors_list2.append([])
            self.neighbors_list3.append([])
        for x in neighbors:
            self.neighbors_list[x[0]].append(x[1])
            self.neighbors_list[x[1]].append(x[0])
        for i in range(self.n_instances):
            for j in self.neighbors_list[i]:
                self.neighbors_list2[i] = list(set(self.neighbors_list2[i]+self.neighbors_list[j]))
        for i in range(self.n_instances):
            for j in self.neighbors_list2[i]:
                self.neighbors_list3[i] = list(set(self.neighbors_list3[i]+self.neighbors_list2[j]))

    def get_distance_from_local_leader(self, neighbors_death):
        leaders = []
        for i in range(self.n_instances):
            if neighbors_death[i] == 0:
                leaders.append(i)
        closest_leader = []
        distance_from_leader = []
        time_from_leader = []
        for i in range(self.n_instances):
            min_distance = math.sqrt(((self.XY[i][0]-self.XY[0][0])**2)+((self.XY[i][1]-self.XY[0][1])**2))
            min_point = 0
            for j in range(len(leaders)):
                new_distance = math.sqrt(((self.XY[i][0]-self.XY[j][0])**2)+((self.XY[i][1]-self.XY[j][1])**2))
                if new_distance < min_distance:
                    min_distance = new_distance
                    min_point = j
            closest_leader.append(min_point)
            distance_from_leader.append(min_distance)
            time_from_leader.append(self.die_times[i] - self.die_times[min_point])

    def identify_leader(self):
        who_is_leader = [-1]*self.n_instances
        for i in range(self.n_instances):
            temp = []
            for k in self.neighbors_list[i]:
                temp.append(who_is_leader[k]) if who_is_leader[k]>-1 else None
            who_is_leader[i] = i if len(temp) == 0 else max(set(temp), key=temp.count)
        neighbors_distance = [-1]*self.n_instances
        for i in range(self.n_instances):
            neighbors_distance[i] = 0 if i == who_is_leader[i] else -1
        for j in range(self.n_instances):
            stop_sign = 0
            for i in range(self.n_instances):
                if neighbors_distance[i] == -1:
                    for k in self.neighbors_list[i]:
                        if neighbors_distance[k] > -1:
                            neighbors_distance[i] = neighbors_distance[k]+1
                            break
                        else:
                            stop_sign += 1
            if stop_sign == 0:
                break
        self.who_is_leader = who_is_leader
        self.death_wave = neighbors_distance

    def get_distance_time_from_leader(self):
        for i in range(self.n_instances):
            time = self.die_times[i] - self.die_times[0]
            distance = math.sqrt(((self.XY[i][0]-self.XY[0][0])**2)+((self.XY[i][1]-self.XY[0][1])**2))
            a = (time, distance)
            self.difference_from_leader.append(a)

    def get_neighbors_stats(self, neighbors_list):
        neighbors_stats_temp = []
        for i in range(self.n_instances):
            neighbors_stats_temp.append([])
            for j in neighbors_list[i]:
                time = self.die_times[i] - self.die_times[j]
                distance = math.sqrt(((self.XY[i][0]-self.XY[j][0])**2)+((self.XY[i][1]-self.XY[j][1])**2))
                a = (j, time, distance)
                neighbors_stats_temp[i].append(a)
        neighbors_stats = []
        for local_list in neighbors_stats_temp:
            neighbors_stats.append(sorted(local_list, key=lambda x: x[2]))
        return neighbors_stats

    def assess_mean(self):
        better_mean = 0
        time_death_means = []
        real_mean_time_death = self.calc_mean(self.neighbors_difference_death_times[0])
        self.mean_time_death = real_mean_time_death * self.time_frame
        for i in range(self.n_scramble):
            temp_mean_time_death = self.calc_mean(self.neighbors_difference_death_times[i + 1])
            time_death_means.append(temp_mean_time_death)
            if temp_mean_time_death > real_mean_time_death:
                better_mean += 1
        time_death_means.sort()
        self.scramble_signficance_95 = time_death_means[int(self.n_scramble * 5 / 100)] * self.time_frame
        self.scramble_signficance_98 = time_death_means[int(self.n_scramble * 2 / 100)] * self.time_frame
        self.scrample_mean_time_death = (sum(time_death_means) / len(time_death_means)) * self.time_frame
        return better_mean / self.n_scramble

    def calc_mean(self, l):
        return sum(l) / float(len(l))

    def draw(self):
        bins = np.linspace(0, 20, 20)
        histogram_scramble = []
        for z in range(len(self.num_blobs)):
            histogram_scramble.append(0)
        for z in range(len(self.neighbors_difference_death_times[0])):
            for i in range(self.n_scramble):
                j = i + 1
                histogram_scramble[self.neighbors_difference_death_times[j][z]]+=1
        avg_scramble_time = []
        for z in range(len(histogram_scramble)):
            histogram_scramble[z] /= self.n_scramble
            n=round(histogram_scramble[z])
            for i in range(n):
                avg_scramble_time.append(z)
        plt.hist(self.neighbors_difference_death_times[0], bins, alpha=0.5, label='Origin')
        plt.hist(avg_scramble_time, bins, alpha=0.5, label="Scramble")
        plt.legend(loc='upper right')
        location = self.pic_path + "\histogram_{}.png".format(self.file_name)
        plt.savefig(location)
        plt.clf()
        plt.plot(self.death_perc_by_time)
        plt.ylabel('Death Perc')
        plt.xlabel('Time of Death')
        location = self.pic_path + "\death_time_{}.png".format(self.file_name)
        plt.savefig(location)
        plt.clf()
        self.draw_time_scatter()
        self.draw_wave_scatter()

    def calc_death_perc(self):
        death_by_time = [0] * len(self.num_blobs)
        for i in range(self.n_instances):
            death_by_time[self.die_times[i]] += 1
        death_by_time_cum = np.cumsum(death_by_time)
        return [x / self.n_instances for x in death_by_time_cum]

    def draw_time_scatter(self):
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        cmap = plt.cm.jet
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
        bounds = np.linspace(0, self.experiment_time, self.experiment_time + 1)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        scat = ax.scatter(self.lisX, self.lisY, c=self.die_times, cmap=cmap, norm=norm)
        cb = plt.colorbar(scat, spacing='proportional', ticks=bounds)
        cb.set_label('Custom cbar')
        ax.set_title('Death Time Scatter')
        location = self.pic_path + "\death_time_scatter_{}.png".format(self.file_name)
        plt.savefig(location)
        plt.clf()

    def draw_wave_scatter(self):
        n=max(self.death_wave)
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        cmap = plt.cm.jet
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
        bounds = np.linspace(0, n, n + 1)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        scat = ax.scatter(self.lisX, self.lisY, c=self.death_wave, cmap=cmap, norm=norm)
        cb = plt.colorbar(scat, spacing='proportional', ticks=bounds)
        cb.set_label('Custom cbar')
        ax.set_title('Death Wave Scatter')
        location = self.pic_path + "\death_wave_scatter_{}.png".format(self.file_name)
        plt.savefig(location)
        plt.clf()


def plot_death_perc_by_time (experiment, file_path):
    for x in experiment:
        plt.plot(x.death_perc_by_time, label=x.file_name)
    plt.ylabel('Death Perc')
    plt.xlabel('Time of Death')
    location = file_path + "Results\All\death_time_all.png"
    plt.savefig(location)
    plt.clf()


def create_combined_file(experiments_list, output_path):
    names = []
    mean_time_death = []
    scramble_signficance_95 = []
    scramble_signficance_98 = []
    scram_mean_time_death = []
    p_value = []
    cell_line = []
    treatment = []
    signal_analyzed = []
    perm_prob = []
    spatial_propagation_index = []
    first = True
    combined_data = None
    for exp in experiments_list:
        if first:
            combined_data = exp.dataframe
            first = False
        else:
            combined_data_temp = combined_data
            combined_data = pd.concat([combined_data_temp, exp.dataframe], ignore_index=True)
        names.append(exp.file_name)
        mean_time_death.append(exp.mean_time_death)
        scramble_signficance_95.append(exp.scramble_signficance_95)
        scramble_signficance_98.append(exp.scramble_signficance_98)
        scram_mean_time_death.append(exp.scrample_mean_time_death)
        p_value.append(1-exp.statistic_score)
        cell_line.append(exp.cell_line)
        treatment.append(exp.treatment)
        signal_analyzed.append(exp.signal_analyzed)
        perm_prob.append(exp.statistic_score)
        propagation_index = ((exp.scramble_signficance_95 - exp.mean_time_death) / exp.scramble_signficance_95)
        spatial_propagation_index.append(propagation_index)

    data = pd.DataFrame(
        {'name': names,
         'mean_time_death': mean_time_death,
         'scramble_signficance ': scramble_signficance_95,
         'scramble_mean_time_death ': scram_mean_time_death,
         'perm_prob ': perm_prob,
         'spatial_propagation_index': spatial_propagation_index,
         'cell_line': cell_line,
         'treatment': treatment,
         'signal_analyzed': signal_analyzed
         })
    location = output_path + "Results/All/"
    if not os.path.exists(location):
        os.makedirs(location)
    data.to_csv(location + "data_combined.csv")
    combined_data.to_csv(location + "data_features_combined.csv")


# single experiment run - exp_file_path should be in csv format
def analyze_single_experiment(exp_file_path, output_path, time_frame, cell_line, treatment, signal_analyzed):
    analysis_results = Analysis(exp_file_path, output_path, time_frame, cell_line, treatment, signal_analyzed)
    # create a combined file with the outputs
    create_combined_file([analysis_results], output_path)


# batch run -  parameters_file_path should be in csv format, exp_folder_path should contain csv files
def analyze_several_experiments(exp_folder_path, output_path, parameters_file_path):
    experiments_data = pd.read_csv(parameters_file_path)  # parameters_file_path should be in csv format
    experiments_list = []
    for file in os.listdir(exp_folder_path):
        # extract the parameters matching to the specific experiment by file name
        file_data = experiments_data[experiments_data['File Name'] == file]
        e_time_frame = int(file_data['Time Interval (min)'].to_numpy())
        e_cell_line = file_data['Cell Line'].to_numpy()[0]
        e_treatment = file_data['Treatment'].to_numpy()[0]
        e_signal_analyzed = file_data['Signal Analyzed'].to_numpy()[0]
        analysis_results = Analysis(exp_folder_path+file, output_path, e_time_frame, e_cell_line, e_treatment, e_signal_analyzed)
        experiments_list.append(analysis_results)
    # create a combined file for the outputs of all experiments
    create_combined_file(experiments_list, output_path)



