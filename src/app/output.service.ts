import { Injectable } from '@angular/core';
import { HttpClient } from "@angular/common/http";
import { Observable } from 'rxjs';
import { Output } from './output'

@Injectable({ providedIn: 'root' })

export class OutputService {

  constructor(private http: HttpClient) { }

  getOutputs(): Observable<Output> {
    return this.http.get<Output>('http://127.0.0.1:5000/');
  }

}