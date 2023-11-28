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

  resetOutputs(): void {
    this.outputs$.next([]);
  }
  
  executeBlock(id: string): void {
    const output: Output = {text: '', other: ''};
    this.http.get<Output>('http://127.0.0.1:5000/'+id+'/').subscribe((o: Output) => {
      output.text = o.text;
      output.other = o.other;
    });
    const outputs = this.outputs$.getValue();
    outputs.push(output);
    this.outputs$.next(outputs);
  }
}
