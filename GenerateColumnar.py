import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
from shapely.geometry import Polygon, box
import random
import math
import time
import os
from dataclasses import dataclass, field
from typing import Literal, Dict

# ==========================================
# 1. Configuration Class
# ==========================================
@dataclass
class GrowthConfig:
    # --- Geometric Domain ---
    domain_width: float 
    domain_height: float 
    thickness_L: float     # Growth thickness (Z-axis)
    
    # --- Mesh Generation ---
    n_elements: int        # Number of grains
    relax_iterations: int  # Mesh regularization steps
    
    # --- Selection Settings ---
    n_joint: int           # Selection denominator (Select 1/n_joint grains)
    
    # --- Tapering (Cone) Settings ---
    enable_cone: bool      # Master switch for tapering
    
    # Cone Angle Mode: 
    # 'fixed'    -> Constant angle for all selected grains
    # 'normal'   -> Gaussian distribution (Mean, Std)
    # 'weighted' -> Discrete probability (Legacy/Specific logic)
    cone_mode: Literal['fixed', 'normal', 'weighted']
    
    # Parameters for the chosen mode
    cone_params: Dict[str, float]

# ==========================================
# 2. Main Simulation Class
# ==========================================
class VoronoiColumnarGrowth:
    def __init__(self, config: GrowthConfig):
        self.cfg = config
        
        # Data Containers
        self.nodes_2d = None
        self.elements = None
        self.element_coords_list = []
        
        # 3D Node Containers
        self.nodes_bottom = None
        self.nodes_top = None
        
        # Selection Records
        self.selected_indices_layer1 = []
        self.selected_indices_layer2 = []

    def _generate_mesh(self):
        """
        Generates Voronoi mesh using Mirroring to ensure perfect boundary filling.
        """
        print(f"[Info] Generating Mesh: N={self.cfg.n_elements}, Domain={self.cfg.domain_width}x{self.cfg.domain_height}")
        w, h = self.cfg.domain_width, self.cfg.domain_height
        domain_box = box(0, 0, w, h)
        
        # Random seeds
        points = np.random.rand(self.cfg.n_elements, 2) * [w, h]
        
        start_time = time.time()
        
        # Lloyd's Relaxation
        for it in range(self.cfg.relax_iterations):
            if (it + 1) % 10 == 0:
                print(f"\r    Relaxation Step: {it+1}/{self.cfg.relax_iterations}", end="")

            # Mirroring (Center, Left, Right, Down, Up)
            pts_c = points
            pts_l = np.copy(pts_c); pts_l[:, 0] = -pts_l[:, 0]
            pts_r = np.copy(pts_c); pts_r[:, 0] = 2*w - pts_r[:, 0]
            pts_d = np.copy(pts_c); pts_d[:, 1] = -pts_d[:, 1]
            pts_u = np.copy(pts_c); pts_u[:, 1] = 2*h - pts_u[:, 1]
            
            points_all = np.vstack([pts_c, pts_l, pts_r, pts_d, pts_u])
            
            try:
                vor = Voronoi(points_all)
            except Exception:
                points += (np.random.rand(*points.shape) - 0.5) * 0.01
                continue

            new_points = []
            for i in range(self.cfg.n_elements):
                region_idx = vor.point_region[i]
                region = vor.regions[region_idx]
                
                if -1 in region or len(region) == 0:
                    new_points.append(points[i])
                    continue
                
                poly_pts = vor.vertices[region]
                poly = Polygon(poly_pts).intersection(domain_box)
                
                if not poly.is_empty and poly.area > 0:
                    new_points.append([poly.centroid.x, poly.centroid.y])
                else:
                    new_points.append(points[i])
            
            points = np.array(new_points)

        print("\n[Info] Relaxation complete. Building topology...")
        
        # Final Generation
        pts_c = points
        pts_l = np.copy(pts_c); pts_l[:, 0] = -pts_l[:, 0]
        pts_r = np.copy(pts_c); pts_r[:, 0] = 2*w - pts_r[:, 0]
        pts_d = np.copy(pts_c); pts_d[:, 1] = -pts_d[:, 1]
        pts_u = np.copy(pts_c); pts_u[:, 1] = 2*h - pts_u[:, 1]
        points_all = np.vstack([pts_c, pts_l, pts_r, pts_d, pts_u])
        
        vor = Voronoi(points_all)
        
        nodes_list = []
        elements_list = []
        element_coords_list = []
        node_map = {}
        next_id = 0
        
        for i in range(self.cfg.n_elements):
            region_idx = vor.point_region[i]
            region = vor.regions[region_idx]
            poly_pts = vor.vertices[region]
            poly = Polygon(poly_pts).intersection(domain_box)
            
            if poly.is_empty: continue
            
            coords = list(poly.exterior.coords)[:-1]
            if poly.area < 0: coords = coords[::-1]
            
            elem_indices = []
            elem_coords_np = []
            
            for x, y in coords:
                key = (round(x, 6), round(y, 6))
                if key not in node_map:
                    node_map[key] = next_id
                    nodes_list.append([x, y])
                    next_id += 1
                elem_indices.append(node_map[key])
                elem_coords_np.append([x, y])
            
            elements_list.append(elem_indices)
            element_coords_list.append(np.array(elem_coords_np))
            
        self.nodes_2d = np.array(nodes_list)
        self.elements = elements_list
        self.element_coords_list = element_coords_list
        print(f"[Info] Mesh constructed with {len(self.elements)} elements.")

    def _get_growth_angle(self):
        """
        Calculates the growth angle based on configuration.
        """
        if not self.cfg.enable_cone:
            return 0.0
            
        mode = self.cfg.cone_mode
        params = self.cfg.cone_params
        angle = 0.0
        
        if mode == 'fixed':
            angle = params.get('value', 0.0)
            
        elif mode == 'normal':
            mu = params.get('mean', 1.0)
            sigma = params.get('std', 0.1)
            angle = np.random.normal(mu, sigma)
            
        elif mode == 'weighted':
            angles = [0.5, 0.9, 1.1, 1.25]
            probs  = [0.45, 0.25, 0.20, 0.10]
            angle = np.random.choice(angles, p=probs)
        
        # Clip to safe positive range
        return max(0.01, min(angle, 15.0))

    def _apply_shrinkage(self, node_array, element_idx, z_height):
        coords = self.element_coords_list[element_idx]
        node_indices = self.elements[element_idx]
        
        if len(coords) == 0: return
        
        centroid = np.mean(coords, axis=0)
        angle = self._get_growth_angle()
        
        # p = L * tan(angle)
        p_val = self.cfg.thickness_L * math.tan(math.radians(angle))
        
        for k, n_id in enumerate(node_indices):
            pt = coords[k]
            dist = np.linalg.norm(pt - centroid)
            
            scale_factor = 1.0
            if dist > 1e-6:
                scale_factor = 1.0 - p_val / dist
            
            # Prevent geometric inversion
            if scale_factor < 0.01: scale_factor = 0.01

            vec = pt - centroid
            new_x = centroid[0] + scale_factor * vec[0]
            new_y = centroid[1] + scale_factor * vec[1]
            
            node_array[n_id, 0] = new_x
            node_array[n_id, 1] = new_y
            node_array[n_id, 2] = z_height

    def process_logic(self):
        # 1. Generate Mesh
        self._generate_mesh()
        
        n_nodes = len(self.nodes_2d)
        self.nodes_bottom = np.zeros((n_nodes, 3))
        self.nodes_top = np.zeros((n_nodes, 3))
        
        # Initialize (Vertical)
        self.nodes_bottom[:, :2] = self.nodes_2d
        self.nodes_bottom[:, 2] = 0.0
        self.nodes_top[:, :2] = self.nodes_2d
        self.nodes_top[:, 2] = self.cfg.thickness_L
        
        # 2. Strict Boundary Detection
        # Any element touching the border is marked as boundary.
        w, h = self.cfg.domain_width, self.cfg.domain_height
        tol = 0.001 
        is_boundary = np.zeros(len(self.elements), dtype=bool)
        
        for i, coords in enumerate(self.element_coords_list):
            if (np.any(coords[:, 0] <= tol) or np.any(coords[:, 0] >= w - tol) or
                np.any(coords[:, 1] <= tol) or np.any(coords[:, 1] >= h - tol)):
                is_boundary[i] = True
        
        # Eligible pool: ONLY Internal grains
        eligible_indices = np.where(~is_boundary)[0].tolist()
        print(f"[Info] Boundary Grains (Fixed): {np.sum(is_boundary)}")
        print(f"[Info] Internal Grains (Eligible for Scaling): {len(eligible_indices)}")
        
        # 3. Layer 1 Selection (Bottom)
        num_select = math.ceil(len(eligible_indices) / self.cfg.n_joint)
        selected_1 = []
        pool = eligible_indices.copy()
        selected_nodes_flat = set()
        
        random.shuffle(pool)
        for idx in pool:
            if len(selected_1) >= num_select: break
            current_nodes = set(self.elements[idx])
            # Avoid immediate neighbors
            if current_nodes.isdisjoint(selected_nodes_flat):
                selected_1.append(idx)
                selected_nodes_flat.update(current_nodes)
        
        self.selected_indices_layer1 = selected_1
        print(f"[Info] Layer 1 Selected (Bottom): {len(selected_1)}")
        
        # Apply Shrinkage Layer 1
        for idx in selected_1:
            self._apply_shrinkage(self.nodes_bottom, idx, 0.0)
            
        # 4. Layer 2 Selection (Top)
        node_to_elems = {}
        for eid, nodes in enumerate(self.elements):
            for nid in nodes:
                if nid not in node_to_elems: node_to_elems[nid] = []
                node_to_elems[nid].append(eid)
        
        adj_map = {i: set() for i in range(len(self.elements))}
        for nid, eids in node_to_elems.items():
            for e1 in eids:
                for e2 in eids:
                    if e1 != e2: adj_map[e1].add(e2)
        
        excluded = set(selected_1)
        for idx in selected_1:
            excluded.update(adj_map[idx])
            
        pool_2 = [i for i in eligible_indices if i not in excluded]
        selected_2 = []
        selected_nodes_flat_2 = set()
        
        random.shuffle(pool_2)
        for idx in pool_2:
            if len(selected_2) >= num_select: break
            current_nodes = set(self.elements[idx])
            if current_nodes.isdisjoint(selected_nodes_flat_2):
                selected_2.append(idx)
                selected_nodes_flat_2.update(current_nodes)
                
        self.selected_indices_layer2 = selected_2
        print(f"[Info] Layer 2 Selected (Top): {len(selected_2)}")
        
        # Apply Shrinkage Layer 2
        for idx in selected_2:
            self._apply_shrinkage(self.nodes_top, idx, self.cfg.thickness_L)

    def export_txt(self, filename):
        current_path = os.getcwd()
        full_path = os.path.join(current_path, filename)
        print(f"[Info] Writing file to: {full_path}")
        
        count = 0
        with open(full_path, 'w') as f:
            for i in range(len(self.elements)):
                elem_indices = self.elements[i]
                n = len(elem_indices)
                for j in range(n):
                    idx_curr = elem_indices[j]
                    idx_next = elem_indices[(j + 1) % n]
                    
                    p1 = self.nodes_bottom[idx_curr]
                    p2 = self.nodes_bottom[idx_next]
                    p3 = self.nodes_top[idx_next]
                    p4 = self.nodes_top[idx_curr]
                    
                    f.write(f"{p1[0]:.8f} {p1[1]:.8f} {p1[2]:.8f}\r\n")
                    f.write(f"{p2[0]:.8f} {p2[1]:.8f} {p2[2]:.8f}\r\n")
                    f.write(f"{p3[0]:.8f} {p3[1]:.8f} {p3[2]:.8f}\r\n")
                    f.write(f"{p4[0]:.8f} {p4[1]:.8f} {p4[2]:.8f}\r\n")
                    count += 1
        print(f"[Success] Written {count} facets.")

    def plot(self):
        print("[Info] Plotting results...")
        plt.figure(figsize=(10, 10))
        c1 = np.array([211, 102, 7]) / 255.0  # Orange
        c2 = np.array([1, 117, 109]) / 255.0  # Cyan
        
        # Plot Bottom
        for i, indices in enumerate(self.elements):
            pts = self.nodes_bottom[indices]
            pts = np.vstack([pts, pts[0]])
            if i in self.selected_indices_layer1:
                plt.fill(pts[:, 0], pts[:, 1], color=c2, alpha=0.8)
            plt.plot(pts[:, 0], pts[:, 1], color=c2, lw=1)
            
        # Plot Top
        for i, indices in enumerate(self.elements):
            pts = self.nodes_top[indices]
            pts = np.vstack([pts, pts[0]])
            if i in self.selected_indices_layer2:
                plt.fill(pts[:, 0], pts[:, 1], color=c1, alpha=0.8)
            plt.plot(pts[:, 0], pts[:, 1], color=c1, lw=1)

        plt.title(f"Simulation: {self.cfg.cone_mode} mode")
        plt.axis('equal')
        plt.xlim(0, self.cfg.domain_width)
        plt.ylim(0, self.cfg.domain_height)
        plt.show()

# ==========================================
# 3. USER CONFIGURATION SECTION
# ==========================================
if __name__ == "__main__":
    
    # --- Define Available Configurations ---
    
    # 1. FIXED ANGLE 
    # Sets a constant 2.0 degrees for all grains
    config_fixed = {
        'enable_cone': True,
        'cone_mode': 'fixed',
        'cone_params': {'value': 2.0}
    }

    # 2. NORMAL DISTRIBUTION
    # Mean = 1.5 deg, Std Dev = 0.5 deg
    config_normal = {
        'enable_cone': True,
        'cone_mode': 'normal',
        'cone_params': {'mean': 1.5, 'std': 0.5}
    }

    # 3. WEIGHTED PROBABILITY 
    # Uses [0.5, 0.9, 1.1, 1.25] with probabilities
    config_weighted = {
        'enable_cone': True,
        'cone_mode': 'weighted',
        'cone_params': {} 
    }
    
    # ---------------------------------------------------------
    # [IMPORTANT] Select your configuration here:
    # ---------------------------------------------------------
    # Change current_settings to config_fixed, config_normal, or config_weighted
    current_settings = config_fixed
    
    # --- Final Config Object ---
    config = GrowthConfig(
        domain_width = 200.0,
        domain_height = 200.0,
        thickness_L = 50.2,
        
        n_elements = 150,      
        relax_iterations = 50,
        
        n_joint = 6, # Determines selection density (1/6 selected)
        
        # Apply the chosen settings
        enable_cone = current_settings['enable_cone'],
        cone_mode   = current_settings['cone_mode'],
        cone_params = current_settings['cone_params']
    )
    
    # --- Run ---
    sim = VoronoiColumnarGrowth(config)
    sim.process_logic()
    sim.export_txt("InitialColumnarIce.txt")
    sim.plot()
