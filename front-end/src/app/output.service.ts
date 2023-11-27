import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { BehaviorSubject, Observable } from 'rxjs';
import { Output } from './output';

@Injectable({
  providedIn: 'root',
})

export class OutputService {
  private readonly outputs$: BehaviorSubject<Output[]> = new BehaviorSubject<Output[]> ([]);
  readonly outputs: Observable<Output[]> = this.outputs$.asObservable();

  constructor(private http: HttpClient) { }

  executeBlocks(): void {
    const output: Output = {
      text: '',
      other: ''
    };
    this.http.get<Output>('http://127.0.0.1:5000/loaddata/').subscribe((o: Output) => {
      output.text = o.text;
      output.other = o.other;
    });
    this.outputs$.next([output]);
    // Execute code and display result on output-display
  }

  getOutputs(): Observable<Output> {
    return this.http.get<Output>('http://127.0.0.1:5000/');
  }
}
