#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Title:       Gravitational Potential of a Homogeneous Tetrahedron
File:        GP_main.py
Description:
    This program computes the gravitational potential of a homogeneous
    tetrahedron at arbitrary observation points (inside, on, and outside)
    using the singularity-free analytical formulation introduced in:
    
        Periyandy, T., & Bevis, M. (2025).
        "The Gravitational Potential Inside, On and Outside of a Homogeneous Tetrahedron."

    The implementation follows the analytical expressions and pseudocode
    described in the paper (Sections 2–3). It provides high numerical stability
    across all spatial regimes through explicit geometric singularity handling.

Authors:
    Thunendran Periyandy  (corresponding author)
        Division of Geodetic Science, School of Earth Sciences,
        The Ohio State University, USA
        Faculty of Geomatics, Sabaragamuwa University of Sri Lanka
        Email: thunendran@gmail.com

    Michael Bevis
        Division of Geodetic Science, School of Earth Sciences,
        The Ohio State University, USA
        Email: mbevis@gmail.com

Version:     1.0
Date:        November 2025
===============================================================================
"""




import numpy as np
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

class TetrahedronPotentialCalculatorOOP:
    def __init__(self, A, B, C, D, P=None, G=1, sigma=1, eps=0):
        """
        Initialize tetrahedron potential calculator.
        If P is provided, compute() will work for single point.
        """
        self.A, self.B, self.C, self.D = self._preserve_orientation_tetrahedron(A, B, C, D)
        self.G = G
        self.sigma = sigma
        self.eps = eps
        self.vertices = {'A': self.A, 'B': self.B, 'C': self.C, 'D': self.D}
        
        # Store P if single computation is needed
        self.P = P
        if P is not None:
            self.face_data = self._precompute_geometric_data_vectorized(P)

    def _preserve_orientation_tetrahedron(self, A, B, C, D):
        V = np.dot(np.cross(B - A, C - A), D - A)
        if V < 0:
            return A, C, B, D
        return A, B, C, D

    def _fast_face_plane_check(self, P, A, B, C):
        n = np.cross(B - A, C - A)
        norm_n = np.linalg.norm(n)
        if norm_n == self.eps:
            return False
        n_unit = n / norm_n
        signed_dist = np.dot(n_unit, P - A)
        return abs(signed_dist) == self.eps

    def _precompute_geometric_data_vectorized(self, P):
        A, B, C, D = self.A, self.B, self.C, self.D
        faces_def = {
            'ACB': (A, C, B), 'ABD': (A, B, D),
            'BCD': (B, C, D), 'CAD': (C, A, D)
        }
        tetra_centroid = (A + B + C + D) / 4.0
        v_map = {id(v): k for k, v in self.vertices.items()}
        vertex_r = {name: np.linalg.norm(v - P) for name, v in self.vertices.items()}

        # Edge data
        edge_lengths, edge_vectors = {}, {}
        for key in ['AB', 'AC', 'AD', 'BC', 'BD', 'CD']:
            v1, v2 = self.vertices[key[0]], self.vertices[key[1]]
            length = np.linalg.norm(v2 - v1)
            edge_lengths[key] = edge_lengths[key[::-1]] = length
            edge_vectors[key] = v2 - v1
            edge_vectors[key[::-1]] = v1 - v2

        # S vectors
        S = {}
        keys = list(self.vertices.keys())
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                vi, vj = self.vertices[keys[i]], self.vertices[keys[j]]
                key = keys[i] + keys[j]
                S[key] = np.cross(vi - P, vj - P)
                S[key[::-1]] = -S[key]

        face_data = {}
        for name, verts in faces_def.items():
            v1, v2, v3 = verts
            if self._fast_face_plane_check(P, v1, v2, v3):
                continue
            n1, n2, n3 = v_map[id(v1)], v_map[id(v2)], v_map[id(v3)]
            n_raw = np.cross(v2 - v1, v3 - v1)
            if np.dot(n_raw, (v1 + v2 + v3) / 3.0 - tetra_centroid) < 0:
                n_raw *= -1
            norm_val = np.linalg.norm(n_raw)
            n_unit = n_raw / norm_val if norm_val > 0 else np.zeros(3)
            e1_key, e2_key, e3_key = n1 + n2, n2 + n3, n3 + n1
            dets = (np.dot(S[e1_key], n_unit),
                    np.dot(S[e2_key], n_unit),
                    np.dot(S[e3_key], n_unit))
            r = (vertex_r[n1], vertex_r[n2], vertex_r[n3])
            edges = (edge_lengths[e1_key], edge_lengths[e2_key], edge_lengths[e3_key])
            dot1 = np.dot(edge_vectors[n1 + n2], edge_vectors[n3 + n1])
            dot2 = np.dot(edge_vectors[n2 + n3], edge_vectors[n1 + n2])
            dot3 = np.dot(edge_vectors[n3 + n1], edge_vectors[n2 + n3])
            diff_z = np.dot(v1 - P, n_unit)
            numerator = diff_z * norm_val
            face_data[name] = {
                'verts': verts,
                'normal': n_unit,
                'determinants': dets,
                'r': r,
                'edges': edges,
                'diff_z': diff_z,
                'numerator': numerator,
                'dots': (dot1, dot2, dot3)
            }
        return face_data

    def _compute_log_terms(self, r, edges):
        r1, r2, r3 = r
        r12, r23, r31 = edges
        def log_term(r_a, r_b, r_ab):
            if r_ab < 0:
                return 0.0
            numerator = r_a + r_b + r_ab
            denominator = r_a + r_b - r_ab
            if denominator <= 0:
                return 0.0
            return np.log(numerator / denominator) / r_ab
        return (
            log_term(r1, r2, r12),
            log_term(r2, r3, r23),
            log_term(r3, r1, r31)
        )

    def _compute_arctangent_terms(self, dets, r, diff_z, numerator, dots):
        det_12, det_23, det_31 = dets
        r1, r2, r3 = r
        dot_BA_AC, dot_CB_BA, dot_AC_CB = dots
        eps = self.eps
        if r1 > eps:
            denom_1 = -((det_31 * det_12) + diff_z**2 * dot_BA_AC) / r1
            S1 = np.arctan2(numerator, denom_1)
        else:
            S1 = 0.0
        if r2 > eps:
            denom_2 = -((det_12 * det_23) + diff_z**2 * dot_CB_BA) / r2
            S2 = np.arctan2(numerator, denom_2)
        else:
            S2 = 0.0
        if r3 > eps:
            denom_3 = -((det_23 * det_31) + diff_z**2 * dot_AC_CB) / r3
            S3 = np.arctan2(numerator, denom_3)
        else:
            S3 = 0.0
        return S1, S2, S3

    def compute(self, P=None):
        """Compute potential at single point."""
        if P is None:
            P = self.P
        face_data = self._precompute_geometric_data_vectorized(P)
        total_U = 0.0
        for data in face_data.values():
            L12, L23, L31 = self._compute_log_terms(data['r'], data['edges'])
            S1, S2, S3 = self._compute_arctangent_terms(
                data['determinants'], data['r'], data['diff_z'], data['numerator'], data['dots']
            )
            det_12, det_23, det_31 = data['determinants']
            diff_z = data['diff_z']
            term1 = diff_z * (det_12 * L12 + det_23 * L23 + det_31 * L31)
            term2 = (diff_z**2) * (S1 + S2 + S3 - np.sign(diff_z) * np.pi)
            total_U += 0.5 * (term1 - term2)
        return self.G * self.sigma * total_U
    
    def compute_vectorized(self, points, n_workers=None):
        """
        Compute potential for multiple points in parallel.

        Parameters
        ----------
        points : (N,3) numpy array
            Evaluation points.
        n_workers : int, optional
            Number of worker processes to use (default: all CPU cores).

        Returns
        -------
        potentials : (N,) numpy array
            Potential at each point.
        """
        if n_workers is None:
            n_workers = multiprocessing.cpu_count()

        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            results = list(executor.map(self.compute, points))

        return np.array(results)

    '''def compute_vectorized(self, points):
        """
        Compute potential for multiple points.
        points: array (N,3)
        Returns: array (N,)
        """
        potentials = np.zeros(points.shape[0])
        for i, P in enumerate(points):
            potentials[i] = self.compute(P)
        return potentials'''
