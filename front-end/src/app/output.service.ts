import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { BehaviorSubject, Observable, firstValueFrom } from 'rxjs';
import { Output } from './output';
import { BlockId } from './block.service';

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
  
  async executeBlock(id: BlockId): Promise<void> {    
    const outputResponse: Output = await firstValueFrom(this.http.get<Output>('http://127.0.0.1:5000/'+id+'/'));
    const outputs = this.outputs$.getValue();
    outputs.push(outputResponse);
    this.outputs$.next(outputs);
  }
}
